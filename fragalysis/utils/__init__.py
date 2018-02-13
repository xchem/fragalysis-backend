from django.contrib.auth.models import User
from rest_framework.authtoken.models import Token


def get_token(request):
    """
    Get the authentication token for a givne request.
    Should just return an un-authenticated user token if nothing.
    :param request:
    :return:
    """
    user = User.objects.get(username=request.user)
    token, created = Token.objects.get_or_create(user=user)
    return token
