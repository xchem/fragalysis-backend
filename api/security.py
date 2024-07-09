# pylint: skip-file
import logging
import os
import threading
from datetime import datetime, timedelta
from functools import cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from wsgiref.util import FileWrapper

from django.conf import settings
from django.db.models import Q
from django.http import Http404, HttpResponse
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
from ispyb.exception import ISPyBConnectionException, ISPyBNoResultException
from rest_framework import viewsets

from viewer.models import Project

from .prometheus_metrics import PrometheusMetrics
from .remote_ispyb_connector import SSHConnector

logger: logging.Logger = logging.getLogger(__name__)


@cache
class CachedContent:
    """
    A static class managing caches proposals/visits for each user.
    Proposals should be collected when has_expired() returns True.
    Content can be written (when the cache for the user has expired)
    and read using the set/get methods.
    """

    _timers: Dict[str, datetime] = {}
    _content: Dict[str, List[str]] = {}
    _cache_period: timedelta = timedelta(
        minutes=settings.SECURITY_CONNECTOR_CACHE_MINUTES
    )
    _cache_lock: threading.Lock = threading.Lock()

    @staticmethod
    def has_expired(username) -> bool:
        assert username
        with CachedContent._cache_lock:
            has_expired = False
            now = datetime.now()
            if username not in CachedContent._timers:
                # User's not known,
                # initialise an entry that will automatically expire
                CachedContent._timers[username] = now
            if CachedContent._timers[username] <= now:
                has_expired = True
                # Expired, reset the expiry time
                CachedContent._timers[username] = now + CachedContent._cache_period
        if has_expired:
            logger.debug("Content expired for '%s'", username)
        return has_expired

    @staticmethod
    def get_content(username):
        with CachedContent._cache_lock:
            if username not in CachedContent._content:
                CachedContent._content[username] = set()
        content = CachedContent._content[username]
        logger.debug("Got content for '%s': %s", username, content)
        return content

    @staticmethod
    def set_content(username, content) -> None:
        with CachedContent._cache_lock:
            CachedContent._content[username] = content.copy()
            logger.debug("Set content for '%s': %s", username, content)


def get_remote_conn(force_error_display=False) -> Optional[SSHConnector]:
    credentials: Dict[str, Any] = {
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
        "ssh_private_key_filename": settings.SSH_PRIVATE_KEY_FILENAME,
        'remote': True,
    }

    credentials.update(**ssh_credentials)

    # Caution: Credentials may not be set in the environment.
    #          Assume the credentials are invalid if there is no host.
    #          If a host is not defined other properties are useless.
    if not credentials["host"]:
        if logging.DEBUG >= logger.level or force_error_display:
            logger.debug("No ISPyB host - cannot return a connector")
        return None

    # Try to get an SSH connection (aware that it might fail)
    logger.debug("Creating remote connector with credentials: %s", credentials)
    conn: Optional[SSHConnector] = None
    try:
        conn = SSHConnector(**credentials)
    except ISPyBConnectionException:
        # The ISPyB connection failed.
        # Nothing else to do here, metrics are already updated
        pass
    except Exception:
        # Any other exception will be a problem with the SSH tunnel connection
        PrometheusMetrics.failed_tunnel()
        if logging.DEBUG >= logger.level or force_error_display:
            logger.info("credentials=%s", credentials)
            logger.exception("Got the following exception creating Connector...")

    if conn:
        logger.debug("Got remote ISPyB connector")
    else:
        logger.debug("Failed to get a remote ISPyB connector")

    return conn


def get_conn(force_error_display=False) -> Optional[Connector]:
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
        if logging.DEBUG >= logger.level or force_error_display:
            logger.info("No ISPyB host - cannot return a connector")
        return None

    logger.info("Creating connector with credentials: %s", credentials)
    conn: Optional[Connector] = None
    try:
        conn = Connector(**credentials)
    except Exception:
        # Log the exception if DEBUG level or lower/finer?
        # The following will not log if the level is set to INFO for example.
        if logging.DEBUG >= logger.level or force_error_display:
            logger.info("credentials=%s", credentials)
            logger.exception("Got the following exception creating Connector...")
    if conn:
        logger.debug("Got connector")
        PrometheusMetrics.new_ispyb_connection()
    else:
        logger.debug("Did not get a connector")
        PrometheusMetrics.failed_ispyb_connection()

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


class ISPyBSafeQuerySet(viewsets.ReadOnlyModelViewSet):
    """
    This ISpyBSafeQuerySet, which inherits from the DRF viewsets.ReadOnlyModelViewSet,
    is used for all views that need to yield (filter) view objects based on a
    user's proposal membership. This requires the view to define the property
    "filter_permissions" to enable this class to navigate to the view object's Project
    (proposal/visit).

    As the ISpyBSafeQuerySet is based on a ReadOnlyModelViewSet, which only provides
    implementations for list() and retrieve() methods, the user will need to provide
    "mixins" for any additional methods the view needs to support (PATCH, PUT, DELETE).
    """

    def get_queryset(self):
        """
        Restricts the returned records to those that belong to proposals
        the user has access to. Without a user only 'open' proposals are returned.
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
        q_filter = self._get_q_filter(proposal_list)
        return self.queryset.filter(q_filter).distinct()

    def get_open_proposals(self):
        """
        Returns the set of proposals anybody can access.
        These consist of any Projects that are marked "open_to_public"
        and any defined via an environment variable.
        """
        open_proposals = set(
            Project.objects.filter(open_to_public=True).values_list("title", flat=True)
        )
        open_proposals.update(settings.PUBLIC_TAS_LIST)
        # Begin Temporary Test Code (1247)
        # Remove any public proposal that's in the restricted list.
        for tas in settings.TEST_RESTRICTED_TAS_LIST:
            if tas in open_proposals:
                open_proposals.remove(tas)
        # End Temporary Test Code (1247)
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
        if CachedContent.has_expired(user.username):
            PrometheusMetrics.new_proposal_cache_miss()
            if conn := get_configured_connector():
                logger.info("Got a connector for '%s'", user.username)
                self._get_proposals_from_connector(user, conn)
            else:
                logger.warning("Failed to get a connector for '%s'", user.username)
        else:
            PrometheusMetrics.new_proposal_cache_hit()

        # The cache has either been updated, has not changed or is empty.
        # Return what we have for the user. Public (open) proposals
        # will be added to what we return if necessary.
        cached_prop_ids = CachedContent.get_content(user.username)
        logger.info(
            "Returning %s cached Proposals for '%s'",
            len(cached_prop_ids),
            user.username,
        )
        return cached_prop_ids

    def _get_proposals_from_connector(self, user, conn):
        """
        Updates the user's proposal cache with the results of a query
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

        # Display the collected results for the user.
        # These will be cached.
        logger.info(
            "%s proposals from %s records for '%s': %s",
            len(prop_id_set),
            len(rs),
            user.username,
            prop_id_set,
        )
        CachedContent.set_content(user.username, prop_id_set)

    def user_is_member_of_target(self, user, target):
        """
        Returns true if the user has access to any proposal the target belongs to.
        """
        target_proposals = [p.title for p in target.project_id.all()]
        user_proposals = self.get_proposals_for_user(user, restrict_to_membership=True)
        return any(proposal in user_proposals for proposal in target_proposals)

    def user_is_member_of_any_given_proposals(self, user, proposals):
        """
        Returns true if the user has access to any proposal in the given
        proposals list. Only one needs to match for permission to be granted.
        We 'restrict_to_membership' to only consider proposals the user
        has explicit membership.
        """
        user_proposals = self.get_proposals_for_user(user, restrict_to_membership=True)
        return any(proposal in user_proposals for proposal in proposals)

    def get_proposals_for_user(self, user, restrict_to_membership=False):
        """
        Returns a list of proposals that the user has access to.

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
            "ispyb_user=%s restrict_to_membership=%s (DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP=%s)",
            ispyb_user,
            restrict_to_membership,
            settings.DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP,
        )
        if ispyb_user:
            if user.is_authenticated:
                logger.debug("Getting proposals from ISPyB...")
                proposals = self._get_proposals_for_user_from_ispyb(user)
        else:
            logger.debug("Getting proposals from Django...")
            proposals = self._get_proposals_for_user_from_django(user)

        # We have all the proposals where the user has authority.
        # Add open/public proposals?
        if (
            not restrict_to_membership
            or settings.DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP
        ):
            proposals.update(self.get_open_proposals())

        # Begin Temporary Test Code (1247)
        # Add test fixed proposals to the given user?
        if user.username in settings.TEST_RESTRICTED_USERS_LIST:
            logger.warning(
                "Adding test restricted proposals for '%s' (%s)",
                user.username,
                settings.TEST_RESTRICTED_TAS_LIST,
            )
            proposals.update(settings.TEST_RESTRICTED_TAS_LIST)
        # End Temporary Test Code (1247)

        # Return the set() as a list()
        return list(proposals)

    def _get_q_filter(self, proposal_list):
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


class ISPyBSafeStaticFiles:
    def get_queryset(self):
        query = ISPyBSafeQuerySet()
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


class ISPyBSafeStaticFiles2(ISPyBSafeStaticFiles):
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
