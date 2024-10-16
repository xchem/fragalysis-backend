import logging

from rest_framework import permissions
from rest_framework.exceptions import PermissionDenied

from api.security import ISPyBSafeQuerySet

_ISPYB_SAFE_QUERY_SET = ISPyBSafeQuerySet()

logger = logging.getLogger(__name__)


class IsObjectProposalMember(permissions.BasePermission):
    """
    Custom permissions to only allow write-access to objects (changes) by users
    who are members of the object's proposals. This permissions class should be used
    in any view that needs to restrict object modifications to users who are members of
    at least one of the object's proposals. This class can be used for objects that
    either have one proposal or many.

    If the object has no proposals, the user is granted access.
    """

    def has_object_permission(self, request, view, obj):
        # Here we check that the user has access to any proposal the object belongs to.
        # Firstly, users must be authenticated
        if not request.user.is_authenticated:
            return False
        # Protect ourselves from views that do not (oddly)
        # have a property called 'filter_permissions'...
        if not hasattr(view, "filter_permissions"):
            raise AttributeError(
                "The view object must define a 'filter_permissions' property"
            )
        # The object's proposal records (one or many) can be obtained via
        # the view's 'filter_permissions' property. A standard
        # django property reference, e.g. 'target__project'.
        object_proposals = []
        attr_value = getattr(obj, view.filter_permissions)

        try:
            attr_value = getattr(obj, view.filter_permissions)
        except AttributeError as exc:
            # Something's gone wrong trying to lookup the project.
            # Log some 'interesting' contextual information...
            logger.info('request=%r', request)
            logger.info('view=%s', view.__class__.__name__)
            logger.info('view.filter_permissions=%s', view.filter_permissions)
            # Get the object's content and dump it for analysis...
            obj_class_name = obj.__class__.__name__
            msg = f"There is no Project at {view.filter_permissions}"
            logger.error(
                "%s - obj=%s vars(base_start_obj)=%s", msg, obj_class_name, vars(obj)
            )
            raise PermissionDenied(msg) from exc

        if attr_value.__class__.__name__ == "ManyRelatedManager":
            # Potential for many proposals...
            object_proposals = [p.title for p in attr_value.all()]
        else:
            # Only one proposal...
            object_proposals = [attr_value.title]
        if not object_proposals:
            raise PermissionDenied(
                detail="Authority cannot be granted - the object is not a part of any Project"
            )
        # Now we have the proposals the object belongs to
        # has the user been associated (in IPSpyB) with any of them?
        if not _ISPYB_SAFE_QUERY_SET.user_is_member_of_any_given_proposals(
            user=request.user, proposals=object_proposals
        ):
            raise PermissionDenied(
                detail="Your authority to access this object has not been given"
            )
        # User is a member of at least one of the object's proposals...
        return True
