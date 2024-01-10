# pylint: skip-file
import logging
import os
import time
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

USER_LIST_DICT: Dict[str, Any] = {}

connector: str = os.environ.get('SECURITY_CONNECTOR', 'ispyb')

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
        "user": os.environ.get("ISPYB_USER"),
        "pw": os.environ.get("ISPYB_PASSWORD"),
        "host": os.environ.get("ISPYB_HOST"),
        "port": os.environ.get("ISPYB_PORT"),
        "db": "ispyb",
        "conn_inactivity": 360,
    }

    ssh_credentials: Dict[str, Any] = {
        'ssh_host': os.environ.get("SSH_HOST"),
        'ssh_user': os.environ.get("SSH_USER"),
        'ssh_password': os.environ.get("SSH_PASSWORD"),
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
        "user": os.environ.get("ISPYB_USER"),
        "pw": os.environ.get("ISPYB_PASSWORD"),
        "host": os.environ.get("ISPYB_HOST"),
        "port": os.environ.get("ISPYB_PORT"),
        "db": "ispyb",
        "conn_inactivity": 360,
    }
    # Caution: Credentials may not have been set in the environment.
    #          Assume the credentials are invalid if there is no host.
    #          If a host is not defined other properties are useless.
    if not credentials["host"]:
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
    if connector == 'ispyb':
        return get_conn()
    elif connector == 'ssh_ispyb':
        return get_remote_conn()
    return None


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

    def get_open_proposals(self):
        """
        Returns the list of proposals anybody can access.
        They are defined via an environment variable
        and are made available as a list of strings (Project titles)
        """
        return settings.PUBLIC_TAS_LIST

    def get_proposals_for_user_from_django(self, user):
        # Get the list of proposals for the user
        if user.pk is None:
            logger.warning("user.pk is None")
            return []
        else:
            prop_ids = list(
                Project.objects.filter(user_id=user.pk).values_list("title", flat=True)
            )
            logger.debug(
                "Got %s proposals for user %s: %s",
                len(prop_ids),
                user.username,
                prop_ids,
            )
            return prop_ids

    def needs_updating(self, user):
        """Returns true of the data collected for a user is out of date.
        It's simple, we just record the last collected timestamp and consider it
        'out of date' (i.e. more than an hour old).
        """
        update_window = 3600
        if user.username not in USER_LIST_DICT:
            USER_LIST_DICT[user.username] = {"RESULTS": [], "TIMESTAMP": 0}
        current_time = time.time()
        if current_time - USER_LIST_DICT[user.username]["TIMESTAMP"] > update_window:
            USER_LIST_DICT[user.username]["TIMESTAMP"] = current_time
            return True
        return False

    def run_query_with_connector(self, conn, user):
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

    def get_proposals_for_user_from_ispyb(self, user):
        # First check if it's updated in the past 1 hour
        needs_updating = self.needs_updating(user)
        logger.info("user=%s needs_updating=%s", user.username, needs_updating)

        if needs_updating:
            conn: Optional[Union[Connector, SSHConnector]] = get_configured_connector()

            # If there is no connection (ISPyB credentials may be missing)
            # then there's nothing we can do except return an empty list.
            # Otherwise run a query for the user.
            if conn is None:
                logger.warning("Failed to get ISPyB connector")
                return []
            rs = self.run_query_with_connector(conn=conn, user=user)
            logger.debug("Connector query rs=%s", rs)

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
            logger.info(
                "Got %s proposals from %s records for user %s: %s",
                len(prop_id_set),
                len(rs),
                user.username,
                prop_id_set,
            )

            # Cache the result and return the result for the user
            USER_LIST_DICT[user.username]["RESULTS"] = list(prop_id_set)
            return USER_LIST_DICT[user.username]["RESULTS"]
        else:
            # Return the previous query (cached for an hour)
            cached_prop_ids = USER_LIST_DICT[user.username]["RESULTS"]
            logger.info(
                "Got %s cached proposals for user %s: %s",
                len(cached_prop_ids),
                user.username,
                cached_prop_ids,
            )
            return cached_prop_ids

    def get_proposals_for_user(self, user):
        """Returns a list of proposals (public and private) that the user has access to."""
        assert user

        proposals = []
        ispyb_user = os.environ.get("ISPYB_USER")
        logger.debug("ispyb_user=%s", ispyb_user)
        if ispyb_user:
            if user.is_authenticated:
                logger.info("Getting proposals from ISPyB...")
                proposals = self.get_proposals_for_user_from_ispyb(user)
            else:
                logger.info(
                    "No proposals (user %s is not authenticated)", user.username
                )
        else:
            logger.info("Getting proposals from Django...")
            proposals = self.get_proposals_for_user_from_django(user)

        # Now add in Target Access Strings that everyone has access to
        # (unless we already have them)
        for open_proposal in self.get_open_proposals():
            if open_proposal not in proposals:
                proposals.append(open_proposal)

        return proposals

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
