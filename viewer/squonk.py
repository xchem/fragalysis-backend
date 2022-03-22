"""
API for Squonk2
Functions that call Squonk2 APIs from Fragalysis
These may ultimately become a separate Pypi library.

"""

def get_type(auth_token):
    """Get the supported File Types

        Returns
        -------
    """

    pass


def put_file(auth_token,
             file,
             as_filename,
             path,
             force):
    """Upload a file into a Project if it doesn't exist

        Returns
        -------
    """
    pass


def post_job_instance(auth_token,
                      project_id,
                      name,
                      callback_url,
                      callback_context,
                      specification):
    """Creates a new job application instance
        callback_url will be of the form:
         <settings.SQUONK_BACKEND_HOST>/api/job_callback/3db328a6-8eeb-4909-a938-aa220b47bfff

        Returns
        -------
    """
    pass