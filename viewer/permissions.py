from rest_framework import permissions
from rest_framework.exceptions import PermissionDenied

from api.security import ISpyBSafeQuerySet

_ISPYB_SAFE_QUERY_SET = ISpyBSafeQuerySet()


class IsManyProposalMember(permissions.BasePermission):
    """Custom permission to only allow objects to be changed
    by users who are members of the object's proposals."""

    def has_object_permission(self, request, view, obj):
        # Here we check that the user has access to one of the proposal(s)
        # the object belongs to. The object's proposal record IDs (Many)
        # is stored in the reference defined in its 'filter_permissions' field.
        object_proposals = [
            p.title for p in getattr(obj, view.filter_permissions).all()
        ]
        # So - is any of the object's proposals in the user's list of proposals?
        if not _ISPYB_SAFE_QUERY_SET.user_is_member_of_any_given_proposals(
            user=request.user, proposals=object_proposals
        ):
            raise PermissionDenied(
                detail="Your authority to access this object has not been given"
            )
        return True
