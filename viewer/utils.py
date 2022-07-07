"""
utils.py

Collection of technical methods tidied up in one location.
"""

import fnmatch
import os
from django.conf import settings

def clean_filename(filepath):
    """Return the "clean" version of a Django filename without the '_abcdefg_' that is created
    when a file is overwritten.

    Given a path and a file, e.g. './media/sdfs/Mpro-x3351_0A_rtEVbqf.sdf' this
    function will return 'Mpro-x3351_0A.sdf'.

    Args:
        filepath

    Returns:
        cleaned filename
    """
    file_split = os.path.splitext(os.path.basename(filepath))
    if fnmatch.fnmatch(
            file_split[0],
            '*_[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'+
            '[a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9]'):
        cleaned_filename = file_split[0][:-8] + file_split[1]
    else:
        cleaned_filename = os.path.basename(filepath)
    return cleaned_filename


def get_https_host(request):
    """Common enabler code for returning https urls

    This is to map links to HTTPS to avoid Mixed Content warnings from Chrome browsers
    SECURE_PROXY_SSL_HEADER is referenced because it is used in redirecting URLs - if
    it is changed it may affect this code.
    Using relative links will probably also work, but This workaround allows both the
    'download structures' button and the DRF API call to work.
    Note that this link will not work on local
    """
    return settings.SECURE_PROXY_SSL_HEADER[1] + '://' + request.get_host()

