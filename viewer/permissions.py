import logging

from rest_framework import permissions

from api.security import ISpyBSafeQuerySet

_LOGGER: logging.Logger = logging.getLogger(__name__)
_ISPYB_SAFE_QUERY_SET: ISpyBSafeQuerySet = ISpyBSafeQuerySet()


class IsProposalMember(permissions.BasePermission):
    """
    Custom permission to only allow owners of an object to edit it.
    """

    def has_object_permission(self, request, view, obj):
        _LOGGER.info("Checking %s", repr(obj))
        _LOGGER.info("view.filter_permissions=%s", view.filter_permissions)
        proposals = _ISPYB_SAFE_QUERY_SET.get_proposals_for_user(
            user=request.user, restrict_to_membership=True
        )
        _LOGGER.info("proposals=%s", proposals)
        return True
