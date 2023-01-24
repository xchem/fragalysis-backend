import logging
import os
import time

from wsgiref.util import FileWrapper
from django.http import Http404
from django.http import HttpResponse
from django.db.models import Q
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
from ispyb.connector.mysqlsp.main import ISPyBNoResultException
from rest_framework import viewsets
from .remote_ispyb_connector import SSHConnector

from viewer.models import Project

logger = logging.getLogger(__name__)

USER_LIST_DICT = {}

connector = os.environ.get('SECURITY_CONNECTOR', 'ispyb')

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


def get_remote_conn():

    ispyb_credentials = {
        "user": os.environ.get("ISPYB_USER"),
        "pw": os.environ.get("ISPYB_PASSWORD"),
        "host": os.environ.get("ISPYB_HOST"),
        "port": os.environ.get("ISPYB_PORT"),
        "db": "ispyb",
        "conn_inactivity": 360,
    }

    ssh_credentials = {
        'ssh_host': os.environ.get("SSH_HOST"),
        'ssh_user': os.environ.get("SSH_USER"),
        'ssh_password': os.environ.get("SSH_PASSWORD"),
        'remote': True
    }

    ispyb_credentials.update(**ssh_credentials)

    # Caution: Credentials may not be set in the environment.
    #          Assume the credentials are invalid if there is no host.
    #          If a host is not defined other properties are useless.
    if not ispyb_credentials["host"]:
        logger.debug("No ISPyB host - cannot return a connector")
        return None

    # Try to get an SSH connection (aware that it might fail)
    conn = None
    try:
        conn = SSHConnector(**ispyb_credentials)
    except ValueError as cex:
        logger.info("ssh_credentials=%s", ssh_credentials)
        logger.error("Got ValueError exception getting SSH connection (%s)", cex)

    return conn


def get_conn():
    credentials = {
        "user": os.environ.get("ISPYB_USER"),
        "pw": os.environ.get("ISPYB_PASSWORD"),
        "host": os.environ.get("ISPYB_HOST"),
        "port": os.environ.get("ISPYB_PORT"),
        "db": "ispyb",
        "conn_inactivity": 360,
    }

    # Caution: Credentials may not be set in the environment.
    #          Assume the credentials are invalid if there is no password
    #          If a host is not defined other properties are useless.
    if not credentials["host"]:
        return None

    conn = Connector(**credentials)
    return conn


class ISpyBSafeQuerySet(viewsets.ReadOnlyModelViewSet):

    def get_queryset(self):
        """
        Optionally restricts the returned purchases to a given proposals
        """
        # The list of proposals this user can have
        proposal_list = self.get_proposals_for_user(self.request.user)
        # Add in the ones everyone has access to
        # (unless we already have it)
        for open_proposal in self.get_open_proposals():
            if open_proposal not in proposal_list:
                proposal_list.append(open_proposal)

        logger.debug('is_authenticated=%s, proposal_list=%s',
                     self.request.user.is_authenticated, proposal_list)

        # Must have a foreign key to a Project for this filter to work.
        # get_q_filter() returns a Q expression for filtering
        q_filter = self.get_q_filter(proposal_list)
        return self.queryset.filter(q_filter).distinct()

    def get_open_proposals(self):
        """
        Returns the list of proposals anybody can access.
        This function is deprecated, instead we should move to the 'open_to_public'
        field rather than using a built-in list of Projects.
        We still add "OPEN" to the list for legacy testing.
        """
        if os.environ.get("TEST_SECURITY_FLAG", False):
            return ["OPEN", "private_dummy_project"]
        else:
            # A list of well-known (built-in) public Projects (Proposals/Visits)
            return ["OPEN", "lb27156"]

    def get_proposals_for_user_from_django(self, user):
        # Get the list of proposals for the user
        if user.pk is None:
            logger.warning("user.pk is None")
            return []
        else:
            prop_ids = list(
                Project.objects.filter(user_id=user.pk).values_list("title", flat=True)
            )
            logger.debug("Got %s proposals: %s", len(prop_ids), prop_ids)
            return prop_ids

    def needs_updating(self, user):
        global USER_LIST_DICT

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
            rs = []
            if conn.server:
                conn.server.stop()
        return rs

    def get_proposals_for_user_from_ispyb(self, user):
        # First check if it's updated in the past 1 hour
        global USER_LIST_DICT

        needs_updating = self.needs_updating(user)
        logger.debug("user=%s needs_updating=%s", user.username, needs_updating)
        
        if needs_updating:
            conn = None
            if connector == 'ispyb':
                conn = get_conn()
            if connector == 'ssh_ispyb':
                conn = get_remote_conn()

            # If there is no connection (ISPyB credentials may be missing)
            # then there's nothing we can do except return an empty list.
            # Otherwise run a query for the user.
            if conn is None:
                logger.warning("Failed to get ISPyB connector")
                return []
            rs = self.run_query_with_connector(conn=conn, user=user)
            logger.debug("Connector query rs=%s", rs)
            
            # Typically you'll find the following fields in each item
            # in the rs response; -
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
            # and one with the 'proposalNumber' and 'sessionNumber' (visits),
            # e.g. ["12345", "12345-1"].
            #
            # These strings would normally correspond to a title value
            # in a Project record.
            #
            # To maintain backward compatibility we create a list of
            # raw proposals and visits and a duplicate set with the 'proposalCode'
            # prefix (converted to upper-case). The proposalCode is 'aa' but 'AA'
            # will be in the project record.
            #
            # e.g. we eventually get this sort of list: -
            # 
            # ["12345", "12345-1", "LB12345", "LB12345-1"]
            prop_id_set = set()
            for record in rs:
                proposal_str = f'{record["proposalNumber"]}'
                visit_str = f'{proposal_str}-{record["sessionNumber"]}'
                prop_id_set.update([proposal_str, visit_str])
                if record["proposalCode"]:
                    proposalCode = str(record["proposalCode"]).upper()
                    prop_id_set.update([f'{proposalCode}{proposal_str}', f'{proposalCode}{visit_str}'])
            logger.debug("Got %s proposals: %s", len(prop_id_set), prop_id_set)

            # Cache the result and return the result for the user
            USER_LIST_DICT[user.username]["RESULTS"] = list(prop_id_set)
            return USER_LIST_DICT[user.username]["RESULTS"]
        else:
            # Return the previous query (cached)
            cached_prop_ids = USER_LIST_DICT[user.username]["RESULTS"]
            logger.debug("Got %s cached proposals: %s", len(cached_prop_ids), cached_prop_ids)
            return cached_prop_ids

    def get_proposals_for_user(self, user):
        """Returns a list of proposals (public and private) that the user has access to.
        """
        assert user
        
        ispyb_user = os.environ.get("ISPYB_USER")
        logger.debug("ispyb_user=%s", ispyb_user)
        if ispyb_user:
            logger.debug("user.is_authenticated=%s", user.is_authenticated)
            if user.is_authenticated:
                return self.get_proposals_for_user_from_ispyb(user)
            else:
                logger.debug("Got no proposals")
                return []
        else:
            return self.get_proposals_for_user_from_django(user)

    def get_q_filter(self, proposal_list):
        """Returns a Q expression representing a (potentially complex) table filter.
        """
        if self.filter_permissions:
            # Q-filter is based on the filter_permissions string
            # whether the resultant Project title in the proposal list
            # OR where the Project is 'open_to_public'
            return Q(**{self.filter_permissions + "__title__in": proposal_list}) |\
                Q(**{self.filter_permissions + "__open_to_public": True})
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
        try:
            queryset = self.get_queryset()
            filter_dict = {self.field_name + "__endswith": self.input_string}
            object = queryset.get(**filter_dict)
            file_name = os.path.basename(str(getattr(object, self.field_name)))

            if hasattr(self, 'file_format'):
                if self.file_format=='raw':
                    file_field = getattr(object, self.field_name)
                    filepath = file_field.path
                    zip_file = open(filepath, 'rb')
                    response = HttpResponse(FileWrapper(zip_file), content_type='application/zip')
                    response['Content-Disposition'] = 'attachment; filename="%s"' % file_name

            else:
                response = HttpResponse()
                response["Content-Type"] = self.content_type
                response["X-Accel-Redirect"] = self.prefix + file_name
                response["Content-Disposition"] = "attachment;filename=" + file_name

            return response
        except Exception:
            raise Http404
