import os
import time

from rest_framework.response import Response
from django.http import Http404
from django.http import HttpResponse
from ispyb.connector.mysqlsp.main import ISPyBMySQLSPConnector as Connector
from ispyb.connector.mysqlsp.main import ISPyBNoResultException
from rest_framework import viewsets

from viewer.models import Project

USER_LIST_DICT = {}


def get_conn():
    credentials = {
        "user": os.environ["ISPYB_USER"],
        "pw": os.environ["ISPYB_PASSWORD"],
        "host": os.environ["ISPYB_HOST"],
        "port": os.environ["ISPYB_PORT"],
        "db": "ispyb",
        "conn_inactivity": 360,
    }
    conn = Connector(**credentials)
    return conn


class ISpyBSafeQuerySet(viewsets.ReadOnlyModelViewSet):

    def get_queryset(self):
        """
        Optionally restricts the returned purchases to a given propsals
        """
        # The list of proposals this user can have
        proposal_list = self.get_proposals_for_user()
        # Add in the ones everyone has access to
        proposal_list.extend(self.get_open_proposals())
        # Must have a directy foreign key (project_id) for it to work
        filter_dict = self.get_filter_dict(proposal_list)
        return self.queryset.filter(**filter_dict).distinct()

    def retrieve(self, request, *args, **kwargs):
        instance = self.get_object()
        serializer = self.get_serializer(instance)
        if isinstance(serializer.data, list):
            count = len(serializer.data)
        else:
            count = len([serializer.data])
        response = {'results': serializer.data, 'count': count, }
        return Response(response)

    def get_open_proposals(self):
        """
        Returns the list of proposals anybody can access
        :return:
        """
        if os.environ.get("TEST_SECURITY_FLAG", False):
            return ["OPEN", "private_dummy_project"]
        else:
            return ["OPEN"]

    def get_proposals_for_user_from_django(self, user):
        # Get the list of proposals for the user
        if user.pk is None:
            return []
        else:
            return list(
                Project.objects.filter(user_id=user.pk).values_list("title", flat=True)
            )

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

    def get_proposals_for_user_from_ispyb(self, user):
        # First check if it's updated in the past 1 hour
        global USER_LIST_DICT
        if self.needs_updating(user):
            with get_conn() as conn:
                core = conn.core
                try:
                    rs = core.retrieve_sessions_for_person_login(user.username)
                except ISPyBNoResultException:
                    rs = []
            visit_ids = [
                str(x["proposalNumber"]) + "-" + str(x["sessionNumber"]) for x in rs
            ]
            prop_ids = [str(x["proposalNumber"]) for x in rs]
            prop_ids.extend(visit_ids)
            USER_LIST_DICT[user.username]["RESULTS"] = prop_ids
            return prop_ids
        else:
            return USER_LIST_DICT[user.username]["RESULTS"]

    def get_proposals_for_user(self):
        user = self.request.user
        get_from_ispyb = os.environ.get("ISPYB_FLAG", True)
        if get_from_ispyb:
            if user.is_authenticated:
                return self.get_proposals_for_user_from_ispyb(user)
            else:
                return []
        else:
            return self.get_proposals_for_user_from_django(user)

    def get_filter_dict(self, proposal_list):
        return {self.filter_permissions + "__title__in": proposal_list}


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
            response = HttpResponse()
            response["Content-Type"] = self.content_type
            response["X-Accel-Redirect"] = self.prefix + file_name
            response["Content-Disposition"] = "attachment;filename=" + file_name
            return response
        except Exception:
            raise Http404
