from rest_framework import permissions
from rest_framework.exceptions import PermissionDenied

from api.security import ISpyBSafeQuerySet

_ISPYB_SAFE_QUERY_SET = ISpyBSafeQuerySet()


class IsManyProposalMember(permissions.BasePermission):
    """
    Custom permission to only allow objects to be changed by users
    who are members of the object's proposals. This permission class should be used
    in any view that need to restrict object modifications to users who are members of
    at least one of the object's proposals.
    """

    def has_object_permission(self, request, view, obj):
        # Users must be authenticated
        if not request.user.is_authenticated:
            return False
        # Here we check that the user has access to any proposal the object belongs to.
        # The object's proposal records (many in this case) are stored in the reference
        # defined in the view's 'filter_permissions' property, e.g. 'project_id'.
        object_proposals = [
            p.title for p in getattr(obj, view.filter_permissions).all()
        ]
        # Has the user been associated (in IPSpyB) with any of the object's proposals?
        if not _ISPYB_SAFE_QUERY_SET.user_is_member_of_any_given_proposals(
            user=request.user, proposals=object_proposals
        ):
            raise PermissionDenied(
                detail="Your authority to access this object has not been given"
            )
        # User is a member of at least one of the object's proposals...
        return True
