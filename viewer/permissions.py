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
        del view

        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in permissions.SAFE_METHODS:
            return True

        _LOGGER.info("Checking %s", repr(obj))
        proposals = _ISPYB_SAFE_QUERY_SET.get_proposals_for_user(request.user)
        _LOGGER.info("Proposals=%s", proposals)
        return True
