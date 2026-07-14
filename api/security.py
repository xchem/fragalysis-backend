# pylint: skip-file
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional
from wsgiref.util import FileWrapper

import ta_auth_connector
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db.models import Q
from django.http import Http404, HttpResponse
from rest_framework import viewsets

from viewer.models import Project, UserRole

from .utils import deployment_mode_is_production

logger: logging.Logger = logging.getLogger(__name__)


def get_restricted_tas_user_proposal(user) -> set[str]:
    """
    Used for debugging access to restricted TAS projects.
    settings.RESTRICTED_TAS_USERS_LIST is a list of strings that
    contain "<user>:<tas>". We inspect this list, and if our user is in it
    we collect and return them.

    This should always return an empty set() in production.
    """
    assert user

    response = set()

    # We ONLY permit the use of RESTRICTED_TAS_USERS
    # when this is not a production deployment
    if not deployment_mode_is_production() and settings.RESTRICTED_TAS_USERS:
        for item in settings.RESTRICTED_TAS_USERS_LIST:
            item_username, item_tas = item.split(':')
            if item_username == user.username:
                response.add(item_tas)

    return response


def ping_configured_connector() -> bool:
    """Pings the connector. If a connection can be obtained it is immediately closed.
    The ping simply provides a way to check the credentials are valid and
    a connection can be made.
    """
    return ta_auth_connector.get_auth_ping().ping == 'OK'


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
            count = len(prop_ids)
            logger.debug(
                "Got %s proposals for '%s': %s",
                count,
                user.username,
                prop_ids,
            )
        return prop_ids

    def user_is_member_of_target(
        self, user, target, restrict_public_to_membership=True
    ):
        """
        Returns true if the user has access to any proposal the target belongs to.
        """
        user_proposals = self.get_proposals_for_user(
            user, restrict_public_to_membership=restrict_public_to_membership
        )
        is_member = target.project.title in user_proposals
        if not is_member:
            logger.warning(
                "Failed membership check user='%s' target='%s' target_proposals=%s",
                user.username,
                target.pk,
                target.project.title,
            )
        return is_member

    def user_is_member_of_any_given_proposals(
        self, user, proposals, restrict_public_to_membership=True
    ):
        """
        Returns true if the user has access to any proposal in the given
        proposals list. Only one needs to match for permission to be granted.
        We 'restrict_public_to_membership' to only consider proposals the user
        has explicit membership.
        """
        user_proposals = self.get_proposals_for_user(
            user, restrict_public_to_membership=restrict_public_to_membership
        )
        is_member = any(proposal in user_proposals for proposal in proposals)
        if not is_member:
            logger.warning(
                "Failed membership check user='%s' proposals=%s",
                user.username,
                proposals,
            )
        return is_member

    def get_proposals_for_user(self, user, restrict_public_to_membership=False):
        """
        Returns a list of proposals that the user has access to.

        If 'restrict_public_to_membership' is set only those proposals/visits where the user
        is a member of the visit will be returned. Otherwise the 'public'
        proposals/visits will also be returned. Typically 'restrict_public_to_membership' is
        used for uploads/changes - this allows us to implement logic that (say)
        only permits explicit members of public proposals to add/load data for that
        project (restrict_public_to_membership=True), but everyone can 'see' public data
        (restrict_public_to_membership=False).
        """
        assert user

        proposals = set()
        ta_auth_service = settings.TA_AUTH_SERVICE
        if ta_auth_service:
            if user.is_authenticated:
                logger.debug(
                    "Getting proposals from TA authenticator (%s)...", ta_auth_service
                )
                proposals = ta_auth_connector.get_auth_target_access(user.username)
            else:
                logger.debug("User is not authenticated")
        else:
            logger.debug("Getting proposals from django...")
            proposals = self._get_proposals_for_user_from_django(user)
        logger.debug("Got %d proposals", len(proposals))

        # We have all the proposals where the user has authority.
        # Add open/public proposals?
        if (
            not restrict_public_to_membership
            or settings.DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP
        ):
            proposals.update(self.get_open_proposals())

        # Finally, add any restricted TAS proposals the user has access to.
        # It uses an environment variable to arbitrarily add proposals for a given user.
        # This is a debug mechanism and should not be used in production.
        # Added during debug effort for 1491.
        proposals.update(get_restricted_tas_user_proposal(user))

        # Return the set() as a list()
        return list(proposals)

    def _get_q_filter(self, proposal_list):
        """Returns a Q expression representing a (potentially complex) table filter."""
        if self.filter_permissions:
            fp = (
                (self.filter_permissions,)
                if isinstance(self.filter_permissions, str)
                else self.filter_permissions
            )
            # Q-filter is based on the filter_permissions strings
            # whether the resultant Project title in the proposal list
            # OR where the Project is 'open_to_public'

            # Chenged during LHS LHS unification to allow multiple
            # permission strings (given as tuple in the
            # view). Important bit: Q objects should be joined with OR

            qry = Q()
            for k in fp:
                qry |= Q(**{k + "__title__in": proposal_list})
                qry |= Q(**{k + "__open_to_public": True})

            return qry

        else:
            # No filter permission?
            # Assume this QuerySet is used for the Project model.
            # Added during 937 development (Access Control).
            #
            # Q-filter is based on the Project title being in the proposal list
            # OR where the Project is 'open_to_public'
            return Q(title__in=proposal_list) | Q(open_to_public=True)


def user_has_loader_role(user) -> bool:
    """True if `user` holds the UserRole.LOADER_ROLE role.

    A Loader bypasses target-access membership checks: they may load data for
    any proposal (see `check_upload_tas_authorisation`) and, symmetrically,
    poll the status of the tasks they start (see `viewer.views.TaskStatusView`).
    """
    return (
        user.is_authenticated and user.roles.filter(name=UserRole.LOADER_ROLE).exists()
    )


@dataclass(frozen=True)
class UploadTASAuthorisationFailure:
    """Why `check_upload_tas_authorisation` refused an upload.

    Exactly one of the two fields is meaningful: `login_required` when the
    user must authenticate first, otherwise `error_body` - the body the view
    should return with a 403 response.
    """

    login_required: bool = False
    error_body: Optional[Dict[str, Any]] = None


def check_upload_tas_authorisation(
    request, target_access_string
) -> Optional[UploadTASAuthorisationFailure]:
    """Authorise an upload/validate request against `target_access_string`.

    Returns None when the user is authorised to proceed, otherwise an
    `UploadTASAuthorisationFailure` the caller must translate into an HTTP
    response (a login redirect or a 403).

    A user holding the UserRole.LOADER_ROLE role bypasses the target-access
    membership check; the bypass is logged as a warning naming the user and
    the role.
    """
    if not settings.AUTHENTICATE_UPLOAD:
        return None

    if request.user.username == 'asap-service':
        logger.warning(
            'Upload attempted with "%s" service account, trying uploader-supplied user',
            request.user.username,
        )
        if 'django-user' in request.headers.keys():
            try:
                user = get_user_model().objects.get(
                    username=request.headers['django-user']
                )
            except get_user_model().DoesNotExist:
                msg = (
                    f'Upload from "{request.user.username}" '
                    + 'service account but fragalysis user not found'
                )
                logger.error(msg)
                return UploadTASAuthorisationFailure(error_body={'error': msg})
        else:
            msg = (
                f'Upload from "{request.user.username}" service '
                'account but fragalysis user not supplied'
            )
            logger.error(msg)
            return UploadTASAuthorisationFailure(error_body={'error': msg})
    else:
        user = request.user

    if not user.is_authenticated:
        return UploadTASAuthorisationFailure(login_required=True)

    if user_has_loader_role(user):
        logger.warning(
            'User "%s" bypassing target-access authorisation for "%s" '
            'via the "%s" role',
            user.username,
            target_access_string,
            UserRole.LOADER_ROLE,
        )
        return None

    proposals = ISPyBSafeQuerySet().get_proposals_for_user(
        user, restrict_public_to_membership=True
    )
    if target_access_string not in proposals:
        logger.warning(
            '(#1712) User %s does not have access to %s (checked %d proposals)',
            user.username,
            target_access_string,
            len(proposals),
        )
        return UploadTASAuthorisationFailure(
            error_body={
                "target_access_string": [
                    f"You are not authorized to upload data to '{target_access_string}'"
                ]
            }
        )

    return None


class ISPyBSafeStaticFiles:
    def get_queryset(self):
        query = ISPyBSafeQuerySet()
        query.request = self.request
        query.filter_permissions = self.permission_string
        query.queryset = self.model.objects.filter()
        queryset = query.get_queryset()
        return queryset

    def get_response(self):
        logger.debug("+ get_response called with: %s", self.input_string)
        try:
            queryset = self.get_queryset()
            filter_dict = {self.field_name + "__endswith": self.input_string}
            logger.debug("filter_dict: %r", filter_dict)
            # instance = queryset.get(**filter_dict)
            instance = queryset.filter(**filter_dict)[0]
            logger.debug("instance: %r", instance)

            file_name = os.path.basename(str(getattr(instance, self.field_name)))

            logger.debug("instance: %r", instance)
            logger.debug("Path to pass to nginx: %s", self.prefix + file_name)

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
        logger.debug("+ get_response called with: %s", self.input_string)
        # it wasn't working because found two objects with test file name
        # so it doesn't help me here..
        try:
            # file_name = Path('/').joinpath(self.prefix).joinpath(self.input_string)
            # file_name = Path(self.prefix).joinpath(self.input_string)
            file_name = str(Path('/').joinpath(self.prefix).joinpath(self.input_string))
            logger.debug("Path to pass to nginx: %s", file_name)
            response = HttpResponse()
            response["Content-Type"] = self.content_type
            response["X-Accel-Redirect"] = file_name
            # response["Content-Disposition"] = "attachment;filename=" + file_name.name
            response["Content-Disposition"] = "attachment;filename=" + self.input_string

            logger.debug("- Response resolved: %r", response)
            return response
        except Exception as exc:
            logger.error(exc, exc_info=True)
            raise Http404 from exc
