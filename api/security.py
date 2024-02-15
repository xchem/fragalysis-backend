# pylint: skip-file
import logging
import os
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, Optional, Union
from wsgiref.util import FileWrapper

from django.conf import settings
from django.db.models import Q
from django.http import Http404, HttpResponse
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
from ispyb.connector.mysqlsp.main import ISPyBNoResultException
from rest_framework import viewsets

from viewer.models import Project

from .remote_ispyb_connector import SSHConnector

logger: logging.Logger = logging.getLogger(__name__)

# Sets of cached query results, indexed by username.
# The cache uses the key 'RESULTS' and the collection time uses the key 'TIMESTAMP'.
# and the time the cache is expires is in 'EXPIRES_AT'
USER_PROPOSAL_CACHE: Dict[str, Dict[str, Any]] = {}
# Period to cache user lists in seconds (on successful reads from the connector)
USER_PROPOSAL_CACHE_MAX_AGE: timedelta = timedelta(
    minutes=settings.SECURITY_CONNECTOR_CACHE_MINUTES
)
# A short period, used when caching of results fails.
# This ensures a rapid retry on failure.
USER_PROPOSAL_CACHE_RETRY_TIMEOUT: timedelta = timedelta(seconds=7)

# example test:
# from rest_framework.test import APIRequestFactory
#
# from rest_framework.test import force_authenticate
# from viewer.views import TargetView
# from django.contrib.auth.models import User
#
# factory = APIRequestFactory()
# view = TargetView.as_view({'get': 'list'})
# user = User.objects.get(username='uzw12877')
# # Make an authenticated request to the view...
# request = factory.get('/api/targets/')
# force_authenticate(request, user=user)
# response = view(request)


def get_remote_conn() -> Optional[SSHConnector]:
    ispyb_credentials: Dict[str, Any] = {
        "user": settings.ISPYB_USER,
        "pw": settings.ISPYB_PASSWORD,
        "host": settings.ISPYB_HOST,
        "port": settings.ISPYB_PORT,
        "db": "ispyb",
        "conn_inactivity": 360,
    }

    ssh_credentials: Dict[str, Any] = {
        'ssh_host': settings.SSH_HOST,
        'ssh_user': settings.SSH_USER,
        'ssh_password': settings.SSH_PASSWORD,
        'remote': True,
    }

    ispyb_credentials.update(**ssh_credentials)

    # Caution: Credentials may not be set in the environment.
    #          Assume the credentials are invalid if there is no host.
    #          If a host is not defined other properties are useless.
    if not ispyb_credentials["host"]:
        logger.debug("No ISPyB host - cannot return a connector")
        return None

    # Try to get an SSH connection (aware that it might fail)
    conn: Optional[SSHConnector] = None
    try:
        conn = SSHConnector(**ispyb_credentials)
    except Exception:
        # Log the exception if DEBUG level or lower/finer?
        # The following wil not log if the level is set to INFO for example.
        if logging.DEBUG >= logger.level:
            logger.info("ispyb_credentials=%s", ispyb_credentials)
            logger.exception("Got the following exception creating SSHConnector...")

    return conn


def get_conn() -> Optional[Connector]:
    credentials: Dict[str, Any] = {
        "user": settings.ISPYB_USER,
        "pw": settings.ISPYB_PASSWORD,
        "host": settings.ISPYB_HOST,
        "port": settings.ISPYB_PORT,
        "db": "ispyb",
        "conn_inactivity": 360,
    }
    # Caution: Credentials may not have been set in the environment.
    #          Assume the credentials are invalid if there is no host.
    #          If a host is not defined other properties are useless.
    if not credentials["host"]:
        logger.debug("No ISPyB host - cannot return a connector")
        return None

    conn: Optional[Connector] = None
    try:
        conn = Connector(**credentials)
    except Exception:
        # Log the exception if DEBUG level or lower/finer?
        # The following wil not log if the level is set to INFO for example.
        if logging.DEBUG >= logger.level:
            logger.info("credentials=%s", credentials)
            logger.exception("Got the following exception creating Connector...")

    return conn


def get_configured_connector() -> Optional[Union[Connector, SSHConnector]]:
    if settings.SECURITY_CONNECTOR == 'ispyb':
        return get_conn()
    elif settings.SECURITY_CONNECTOR == 'ssh_ispyb':
        return get_remote_conn()
    return None


def ping_configured_connector() -> bool:
    """Pings the connector. If a connection can be obtained it is immediately closed.
    The ping simply provides a way to check the credentials are valid and
    a connection can be made.
    """
    conn: Optional[Union[Connector, SSHConnector]] = None
    if settings.SECURITY_CONNECTOR == 'ispyb':
        conn = get_conn()
    elif settings.SECURITY_CONNECTOR == 'ssh_ispyb':
        conn = get_remote_conn()
        if conn is not None:
            conn.stop()
    return conn is not None


class ISpyBSafeQuerySet(viewsets.ReadOnlyModelViewSet):
    def get_queryset(self):
        """
        Optionally restricts the returned purchases to a given proposals
        """
        # The list of proposals this user can have
        proposal_list = self.get_proposals_for_user(self.request.user)
        logger.debug(
            'is_authenticated=%s, proposal_list=%s',
            self.request.user.is_authenticated,
            proposal_list,
        )

        # Must have a foreign key to a Project for this filter to work.
        # get_q_filter() returns a Q expression for filtering
        q_filter = self.get_q_filter(proposal_list)
        return self.queryset.filter(q_filter).distinct()

    def _get_open_proposals(self):
        """
        Returns the set of proposals anybody can access.
        These consist of any Projects that are marked "open_to_public"
        and any defined via an environment variable.
        """
        open_proposals = set(
            Project.objects.filter(open_to_public=True).values_list("title", flat=True)
        )
        open_proposals.update(settings.PUBLIC_TAS_LIST)
        return open_proposals

    def _get_proposals_for_user_from_django(self, user):
        prop_ids = set()
        # Get the set() of proposals for the user
        if user.pk is None:
            logger.warning("user.pk is None")
        else:
            prop_ids.update(
                Project.objects.filter(user_id=user.pk).values_list("title", flat=True)
            )
            logger.info(
                "Got %s proposals for '%s': %s",
                len(prop_ids),
                user.username,
                prop_ids,
            )
        return prop_ids

    def _cache_needs_updating(self, user):
        """True of the data for a user now needs to be collected
        (e.g. the cache is out of date). The response is also True for the first
        call for each user. When data is successfully collected you need to
        call '_populate_cache()' with the user and new cache content.
        This will set the cache content and the cache timestamp.
        """
        now = datetime.now()
        if user.username not in USER_PROPOSAL_CACHE:
            # Unknown user - initilise the entry for this user.
            # And make suer it immediately expires!
            USER_PROPOSAL_CACHE[user.username] = {
                "RESULTS": set(),
                "TIMESTAMP": None,
                "EXPIRES_AT": now,
            }

        # Has the cache expired?
        return now >= USER_PROPOSAL_CACHE[user.username]["EXPIRES_AT"]

    def _populate_cache(self, user, new_content):
        """Called by code that collects content to replace the cache with new content,
        this is typically from '_get_proposals_from_connector()'. The underlying map's
        TIMESTAMP for the user will also be set (to 'now') to reflect the time the
        cache was most recently populated.
        """
        username = user.username
        USER_PROPOSAL_CACHE[username]["RESULTS"] = new_content.copy()
        # Set the timestamp (which records when the cache was populated with 'stuff')
        # and ensure it will now expire after USER_PROPOSAL_CACHE_SECONDS.
        now = datetime.now()
        USER_PROPOSAL_CACHE[username]["TIMESTAMP"] = now
        USER_PROPOSAL_CACHE[username]["EXPIRES_AT"] = now + USER_PROPOSAL_CACHE_MAX_AGE
        logger.info(
            "USER_PROPOSAL_CACHE populated for '%s' (expires at %s)",
            username,
            USER_PROPOSAL_CACHE[username]["EXPIRES_AT"],
        )

    def _mark_cache_collection_failure(self, user):
        """Called by code that collects content to indicate that although the cache
        should have been collected it has not (trough some other problem).
        Under these circumstances the cache will not be updated but we have the opportunity
        to set a new, short, 'expiry' time. In this way, cache collection will occur
        again soon. The cache and its timestamp are left intact.
        """
        now = datetime.now()
        USER_PROPOSAL_CACHE[user.username]["EXPIRES_AT"] = (
            now + USER_PROPOSAL_CACHE_RETRY_TIMEOUT
        )

    def _run_query_with_connector(self, conn, user):
        core = conn.core
        try:
            rs = core.retrieve_sessions_for_person_login(user.username)
            if conn.server:
                conn.server.stop()
        except ISPyBNoResultException:
            logger.warning("No results for user=%s", user.username)
            rs = []
            if conn.server:
                conn.server.stop()
        return rs

    def _get_proposals_for_user_from_ispyb(self, user):
        if self._cache_needs_updating(user):
            logger.info("user='%s' (needs_updating)", user.username)
            if conn := get_configured_connector():
                logger.debug("Got a connector for '%s'", user.username)
                self._get_proposals_from_connector(user, conn)
            else:
                logger.warning("Failed to get a connector for '%s'", user.username)
                self._mark_cache_collection_failure(user)

        # The cache has either been updated, has not changed or is empty.
        # Return what we have for the user. If required, public (open) proposals
        # will be added to what we return.
        cached_prop_ids = USER_PROPOSAL_CACHE[user.username]["RESULTS"]
        logger.info(
            "Got %s proposals for '%s': %s",
            len(cached_prop_ids),
            user.username,
            cached_prop_ids,
        )
        return cached_prop_ids

    def _get_proposals_from_connector(self, user, conn):
        """Updates the USER_LIST_DICT with the results of a query
        and marks it as populated.
        """
        assert user
        assert conn

        rs = self._run_query_with_connector(conn=conn, user=user)

        # Typically you'll find the following fields in each item
        # in the rs response: -
        #
        #    'id': 0000000,
        #    'proposalId': 00000,
        #    'startDate': datetime.datetime(2022, 12, 1, 15, 56, 30)
        #    'endDate': datetime.datetime(2022, 12, 3, 18, 34, 9)
        #    'beamline': 'i00-0'
        #    'proposalCode': 'lb'
        #    'proposalNumber': '12345'
        #    'sessionNumber': 1
        #    'comments': None
        #    'personRoleOnSession': 'Data Access'
        #    'personRemoteOnSession': 1
        #
        # Iterate through the response and return the 'proposalNumber' (proposals)
        # and one with the 'proposalNumber' and 'sessionNumber' (visits), each
        # prefixed by the `proposalCode` (if present).
        #
        # Codes are expected to consist of 2 letters.
        # Typically: lb, mx, nt, nr, bi
        #
        # These strings should correspond to a title value in a Project record.
        # and should get this sort of list: -
        #
        # ["lb12345", "lb12345-1"]
        #              --      -
        #              | ----- |
        #           Code   |   Session
        #               Proposal
        prop_id_set = set()
        for record in rs:
            pc_str = ""
            if "proposalCode" in record and record["proposalCode"]:
                pc_str = f'{record["proposalCode"]}'
            pn_str = f'{record["proposalNumber"]}'
            sn_str = f'{record["sessionNumber"]}'
            proposal_str = f'{pc_str}{pn_str}'
            proposal_visit_str = f'{proposal_str}-{sn_str}'
            prop_id_set.update([proposal_str, proposal_visit_str])

        # Always display the collected results for the user.
        # These will be cached.
        logger.debug(
            "%s proposals from %s records for '%s': %s",
            len(prop_id_set),
            len(rs),
            user.username,
            prop_id_set,
        )

        # Replace the cache with what we've collected
        self._populate_cache(user, prop_id_set)

    def get_proposals_for_user(self, user, restrict_to_membership=False):
        """Returns a list of proposals that the user has access to.

        If 'restrict_to_membership' is set only those proposals/visits where the user
        is a member of the visit will be returned. Otherwise the 'public'
        proposals/visits will also be returned. Typically 'restrict_to_membership' is
        used for uploads/changes - this allows us to implement logic that (say)
        only permits explicit members of public proposals to add/load data for that
        project (restrict_to_membership=True), but everyone can 'see' public data
        (restrict_to_membership=False).
        """
        assert user

        proposals = set()
        ispyb_user = settings.ISPYB_USER
        logger.debug(
            "ispyb_user=%s restrict_to_membership=%s",
            ispyb_user,
            restrict_to_membership,
        )
        if ispyb_user:
            if user.is_authenticated:
                logger.info("Getting proposals from ISPyB...")
                proposals = self._get_proposals_for_user_from_ispyb(user)
        else:
            logger.info("Getting proposals from Django...")
            proposals = self._get_proposals_for_user_from_django(user)

        # We have all the proposals where the user has authority.
        # Add open/public proposals?
        if not restrict_to_membership:
            proposals.update(self._get_open_proposals())

        # Return the set() as a list()
        return list(proposals)

    def get_q_filter(self, proposal_list):
        """Returns a Q expression representing a (potentially complex) table filter."""
        if self.filter_permissions:
            # Q-filter is based on the filter_permissions string
            # whether the resultant Project title in the proposal list
            # OR where the Project is 'open_to_public'
            return Q(**{self.filter_permissions + "__title__in": proposal_list}) | Q(
                **{self.filter_permissions + "__open_to_public": True}
            )
        else:
            # No filter permission?
            # Assume this QuerySet is used for the Project model.
            # Added during 937 development (Access Control).
            #
            # Q-filter is based on the Project title being in the proposal list
            # OR where the Project is 'open_to_public'
            return Q(title__in=proposal_list) | Q(open_to_public=True)


class ISpyBSafeStaticFiles:
    def get_queryset(self):
        query = ISpyBSafeQuerySet()
        query.request = self.request
        query.filter_permissions = self.permission_string
        query.queryset = self.model.objects.filter()
        queryset = query.get_queryset()
        return queryset

    def get_response(self):
        logger.info("+ get_response called with: %s", self.input_string)
        try:
            queryset = self.get_queryset()
            filter_dict = {self.field_name + "__endswith": self.input_string}
            logger.info("filter_dict: %r", filter_dict)
            # instance = queryset.get(**filter_dict)
            instance = queryset.filter(**filter_dict)[0]
            logger.info("instance: %r", instance)

            file_name = os.path.basename(str(getattr(instance, self.field_name)))

            logger.info("instance: %r", instance)
            logger.info("Path to pass to nginx: %s", self.prefix + file_name)

            if hasattr(self, 'file_format'):
                if self.file_format == 'raw':
                    file_field = getattr(object, self.field_name)
                    filepath = file_field.path
                    zip_file = open(filepath, 'rb')
                    response = HttpResponse(
                        FileWrapper(zip_file), content_type='application/zip'
                    )
                    response['Content-Disposition'] = (
                        'attachment; filename="%s"' % file_name
                    )

            else:
                response = HttpResponse()
                response["Content-Type"] = self.content_type
                response["X-Accel-Redirect"] = self.prefix + file_name
                response["Content-Disposition"] = "attachment;filename=" + file_name

            return response
        except Exception as exc:
            logger.error(exc, exc_info=True)
            raise Http404 from exc


class ISpyBSafeStaticFiles2(ISpyBSafeStaticFiles):
    def get_response(self):
        logger.info("+ get_response called with: %s", self.input_string)
        # it wasn't working because found two objects with test file name
        # so it doesn't help me here..
        try:
            # file_name = Path('/').joinpath(self.prefix).joinpath(self.input_string)
            # file_name = Path(self.prefix).joinpath(self.input_string)
            file_name = str(Path('/').joinpath(self.prefix).joinpath(self.input_string))
            logger.info("Path to pass to nginx: %s", file_name)
            response = HttpResponse()
            response["Content-Type"] = self.content_type
            response["X-Accel-Redirect"] = file_name
            # response["Content-Disposition"] = "attachment;filename=" + file_name.name
            response["Content-Disposition"] = "attachment;filename=" + self.input_string

            logger.info("- Response resolved: %r", response)
            return response
        except Exception as exc:
            logger.error(exc, exc_info=True)
            raise Http404 from exc
