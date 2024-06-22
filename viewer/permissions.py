import logging
from typing import List

from rest_framework import permissions

from api.security import ISpyBSafeQuerySet

_LOGGER: logging.Logger = logging.getLogger(__name__)
_ISPYB_SAFE_QUERY_SET: ISpyBSafeQuerySet = ISpyBSafeQuerySet()


class IsProposalMember(permissions.BasePermission):
    """
    Custom permission to only allow owners of an object to edit it.
    """

    def has_object_permission(self, request, view, obj):
        # Here we check that the user has access to the proposal(s)
        # the object belongs to.
        # What proposals does the Object belong to?
        # The object's proposal record ID is stored in the filter_permissions field.
        object_proposals: List[str] = [
            p.title for p in getattr(obj, view.filter_permissions).all()
        ]
        _LOGGER.info("obj.proposals=%s", object_proposals)
        # What proposals does the user have access to?
        # We restrict the query to only those proposals
        # where the user has explicit membership
        # ensuring that the user is a member of an open proposal in order to access it.
        user_proposals = _ISPYB_SAFE_QUERY_SET.get_proposals_for_user(user=request.user)
        _LOGGER.info("user_proposals=%s", user_proposals)
        # So - is any of the object's proposals in the user's list of proposals?
        # Only one needs to match for permission to be granted.
        permission = any(proposal in user_proposals for proposal in object_proposals)
        _LOGGER.info("permission=%s", permission)
        return permission
