"""An agent for the Squonk2 (Data Manger and Account Server) API.
This module 'simplifies' the use of the Squonk2 Python client package.
"""

# Refer to the accompanying low-level-design document: -
# https://docs.google.com/document/d/1lFpN29dK1luz80lwSGi0Rnj1Rqula2_CTRuyWDUBu14

from collections import namedtuple
import logging
import os
from typing import List, Optional, Tuple
from urllib.parse import ParseResult, urlparse
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings

from django.core.exceptions import ObjectDoesNotExist
from squonk2.auth import Auth
from squonk2.as_api import AsApi, AsApiRv
from squonk2.dm_api import DmApi, DmApiRv
import requests
from requests import Response
from wrapt import synchronized

from api.security import ISpyBSafeQuerySet
from viewer.models import User, SessionProject, Project
from viewer.models import Squonk2Project, Squonk2Org, Squonk2Unit

_LOGGER: logging.Logger = logging.getLogger(__name__)

# Response value for the agent methods
Squonk2AgentRv: namedtuple = namedtuple('Squonk2AgentRv', ['success', 'msg'])
SuccessRv: Squonk2AgentRv = Squonk2AgentRv(success=True, msg=None)

# Parameters common to each named Tuple
CommonParams: namedtuple = namedtuple("CommonParams",
                                      ["user_id",
                                       "access_id",
                                       "target_id",
                                       "session_id"])

# Named tuples are used to pass parameters to the agent methods.
# RunJob, used in run_job()
RunJobParams: namedtuple = namedtuple("RunJob", ["common",
                                                 "job_spec",
                                                 "callback_url"])

# Send, used in send()
# Access ID is the Fragalysis Project record ID
# Session ID is the Fragalysis SessionProject record ID
# Target ID is the Fragalysis Target record ID
# Snapshot ID is the Fragalysis Snapshot record ID
SendParams: namedtuple = namedtuple("Send", ["common",
                                             "snapshot_id"])

_SUPPORTED_PRODUCT_FLAVOURS: List[str] = ["BRONZE", "SILVER", "GOLD"]

# Squonk2 Have defined limits - assumed here.
# verify with your Squonk2 installation.
# How long are Squonk2 'names'?
_SQ2_MAX_NAME_LENGTH: int = 80

# A slug used for names this Fragalysis will create
# and a prefix string. So Squonk2 objects will be called "Fragalysis {slug}"
_MAX_SLUG_LENGTH: int = 10
_SQ2_NAME_PREFIX: str = "Fragalysis"

# Built-in
_SQ2_PRODUCT_TYPE: str = 'DATA_MANAGER_PROJECT_TIER_SUBSCRIPTION'

# True if the code's in Test Mode
_TEST_MODE: bool = False


class Squonk2Agent:
    """Helper class that simplifies access to the Squonk2 Python client.
    Users shouldn't instantiate the class directly, instead they should
    get access to the class singleton via a call to 'get_squonk2_agent()'.

    The class methods protect the caller from using them unless a) the class has
    sufficient configuration and b) the Squonk2 services are 'alive'.
    """

    def __init__(self):
        """Initialise the instance, loading from the environment.
        """

        # Primary configuration of the module is via the container environment.
        # We need to recognise that some or all of these may not be defined.
        # All run-time config that's required is given a __CFG prefix to
        # simplify checking whether all that's required has been defined.
        #
        # The SQUONK2_SLUG is limited to 10 characters, when combined with
        # "Fragalysis {SLUG} ", this leaves (80-22) 58 characters for the
        # use with the target-access-string and session project strings
        # to form Squonk2 Unit and Project names.
        self.__CFG_SQUONK2_ASAPI_URL: Optional[str] =\
            os.environ.get('SQUONK2_ASAPI_URL')
        self.__CFG_SQUONK2_DMAPI_URL: Optional[str] =\
            os.environ.get('SQUONK2_DMAPI_URL')
        self.__CFG_SQUONK2_UI_URL: Optional[str] =\
            os.environ.get('SQUONK2_UI_URL')
        self.__CFG_SQUONK2_ORG_UUID: Optional[str] =\
            os.environ.get('SQUONK2_ORG_UUID')
        self.__CFG_SQUONK2_UNIT_BILLING_DAY: Optional[str] =\
            os.environ.get('SQUONK2_UNIT_BILLING_DAY')
        self.__CFG_SQUONK2_PRODUCT_FLAVOUR: Optional[str] =\
            os.environ.get('SQUONK2_PRODUCT_FLAVOUR')
        self.__CFG_SQUONK2_SLUG: Optional[str] =\
            os.environ.get('SQUONK2_SLUG', '')[:_MAX_SLUG_LENGTH]
        self.__CFG_SQUONK2_ORG_OWNER: Optional[str] =\
            os.environ.get('SQUONK2_ORG_OWNER')
        self.__CFG_SQUONK2_ORG_OWNER_PASSWORD: Optional[str] =\
            os.environ.get('SQUONK2_ORG_OWNER_PASSWORD')
        self.__CFG_OIDC_AS_CLIENT_ID: Optional[str] = \
            os.environ.get('OIDC_AS_CLIENT_ID')
        self.__CFG_OIDC_DM_CLIENT_ID: Optional[str] = \
            os.environ.get('OIDC_DM_CLIENT_ID')
        self.__CFG_OIDC_KEYCLOAK_REALM: Optional[str] = \
            os.environ.get('OIDC_KEYCLOAK_REALM')

        # Optional config (no '__CFG_' prefix)
        self.__DUMMY_SESSION_TITLE: Optional[str] =\
            os.environ.get('DUMMY_SESSION_TITLE')
        self.__DUMMY_USER: Optional[str] =\
            os.environ.get('DUMMY_USER')
        self.__DUMMY_TAS: Optional[str] =\
            os.environ.get('DUMMY_TAS')
        self.__SQUONK2_VERIFY_CERTIFICATES: Optional[str] = \
            os.environ.get('SQUONK2_VERIFY_CERTIFICATES')

        # The integer billing day, valid if greater than zero
        self.__unit_billing_day: int = 0
        # True if configured...
        self.__configuration_checked: bool = False
        self.__configured: bool = False
        # OIDC hostname and realm.
        # Extracted during configuration check from the OIDC variable
        self.__oidc_hostname: str = ''
        self.__oidc_realm: str = ''
        # Ignore cert errors? (no)
        self.__verify_certificates: bool = True

        # Set when pre-flight checks have passed.
        # When they've been done we can safely (?) continue to use the
        # Squonk2 Python client.
        self.__pre_flight_check_status: bool = False
        # The record ID of the Squonk2Org for this deployment.
        # Set on successful 'pre-flight-check'
        self.__org_record: Optional[Squonk2Org] = None

        self.__org_owner_as_token: str = ''
        self.__org_owner_dm_token: str = ''
        self.__keycloak_hostname: str = ''
        self.__keycloak_realm: str = ''

        # The Safe QuerySet from the security module.
        # Used when we are given an access_id.
        # It allows us to check that a user is permitted to use the access ID
        # and relies on ISPyB credentials present in the environment.
        self.__ispyb_safe_query_set: ISpyBSafeQuerySet = ISpyBSafeQuerySet()

    def _get_user_name(self, user_id: int) -> str:
        # Gets the username (if id looks sensible)
        # or a fixed value (used for testing)
        if user_id == 0:
            assert _TEST_MODE
            _LOGGER.warning('Caution - in TEST mode, using __DUMMY_USER (%s)',
                            self.__DUMMY_USER)

        user_name: str = self.__DUMMY_USER
        if user_id:
            user: User  = User.objects.filter(id=user_id).first()
            assert user
            user_name = user.username
        assert user_name
        return user_name

    def _get_target_access_string(self, access_id: int) -> str:
        # Gets the Target Access String (if id looks sensible)
        # or a fixed value (used for testing)
        if access_id == 0:
            assert _TEST_MODE
            _LOGGER.warning('Caution - in TEST mode, using __DUMMY_TAS (%s)',
                            self.__DUMMY_TAS)

        # Get the "target access string" (test mode or otherwise)
        target_access_string: str = self.__DUMMY_TAS
        if access_id:
            project: Optional[Project] = Project.objects.filter(id=access_id).first()
            assert project
            target_access_string = project.title
        assert target_access_string
        return target_access_string

    def _get_session_title(self, session_id: int) -> str:
        # Gets the Session title (if id looks sensible)
        # or a fixed value (used for testing)
        if session_id == 0:
            assert _TEST_MODE
            _LOGGER.warning('Caution - in TEST mode, using __DUMMY_TAS__DUMMY_SESSION_TITLE (%s)',
                            self.__DUMMY_SESSION_TITLE)

        session_title: str = self.__DUMMY_SESSION_TITLE
        if session_id:
            session_project: SessionProject = SessionProject.objects.filter(id=params.session_id).first()
            assert session_project
            session_title = session_project.title
        assert session_title
        return session_title

    def _build_unit_name(self, target_access_string: str) -> Tuple[str, str]:
        assert target_access_string
        # AS Units are named using the Target Access String (TAS)
        # which, in Diamond will be a "visit" string like "lb00000-1"
        name: str = f'{_SQ2_NAME_PREFIX} {self.__CFG_SQUONK2_SLUG} /{target_access_string}/'
        return name[:_SQ2_MAX_NAME_LENGTH], name

    def _build_product_name(self, username: str, session_string: str) -> Tuple[str, str]:
        """Builds a Product name, returning the truncated and un-truncated form"""
        assert username
        assert session_string
        # AS Products (there's a 1:1 mapping to DM Projects)
        # are named using the user and the session

        # The Product name characters are not restricted
        identifier: str = f'{username}::{session_string}'
        name: str = f'{_SQ2_NAME_PREFIX} {self.__CFG_SQUONK2_SLUG} {identifier}'
        return name[:_SQ2_MAX_NAME_LENGTH], name

    def _build_project_name(self, user_id: int, session_id: int) -> Tuple[str, str]:
        assert user_id
        assert session_id
        # DM Projects (there's a 1:1 mapping to Products)
        # are named using the user and the session

        # The Project name characters are RESTRICTED,
        # and need to be limited to characters that are
        # valid for use with RFC 1123 Label Names
        identifier: str = f'{user_id}-{session_id}'
        name: str = f'{_SQ2_NAME_PREFIX} {self.__CFG_SQUONK2_SLUG} {identifier}'
        return name[:_SQ2_MAX_NAME_LENGTH], name

    def _get_squonk2_owner_tokens(self) -> Optional[Tuple[str, str]]:
        """Gets access tokens for the Squonk2 organisation owner.
        This sets the __keycloak_hostname member and also returns the tokens,
        getting one for the AS and one for the DM.
        """
        assert self.__keycloak_hostname

        _LOGGER.debug('__keycloak_hostname="%s" __keycloak_realm="%s"'
                      ' dm-client=%s as-client=%s org_owner=%s',
                     self.__keycloak_hostname,
                     self.__keycloak_realm,
                     self.__CFG_OIDC_DM_CLIENT_ID,
                     self.__CFG_OIDC_AS_CLIENT_ID,
                     self.__CFG_OIDC_AS_CLIENT_ID)

        self.__org_owner_as_token = Auth.get_access_token(
            keycloak_url="https://" + self.__keycloak_hostname + "/auth",
            keycloak_realm=self.__keycloak_realm,
            keycloak_client_id=self.__CFG_OIDC_AS_CLIENT_ID,
            username=self.__CFG_SQUONK2_ORG_OWNER,
            password=self.__CFG_SQUONK2_ORG_OWNER_PASSWORD,
        )
        if not self.__org_owner_as_token:
            _LOGGER.warning('Failed to get access token for AS Organisation owner')
            return None

        self.__org_owner_dm_token = Auth.get_access_token(
            keycloak_url="https://" + self.__keycloak_hostname + "/auth",
            keycloak_realm=self.__keycloak_realm,
            keycloak_client_id=self.__CFG_OIDC_DM_CLIENT_ID,
            username=self.__CFG_SQUONK2_ORG_OWNER,
            password=self.__CFG_SQUONK2_ORG_OWNER_PASSWORD,
        )
        if not self.__org_owner_dm_token:
            _LOGGER.warning('Failed to get access token for DM as AS Organisation owner')
            return None

        # OK if we get here
        return self.__org_owner_as_token, self.__org_owner_dm_token

    def _pre_flight_checks(self) -> Squonk2AgentRv:
        """Execute pre-flight checks,
        can be called multiple times, it acts only once.
        """

        # If a Squonk2Org record exists its UUID cannot have changed.
        # We cannot change the organisation once deployed. The corresponding Units,
        # Products and Projects are organisation-specific. The Squonk2Org table
        # records the organisation ID and the Account Server URL where the ID
        # is valid. None of these values can change once deployed.

        squonk2_org: Optional[Squonk2Org] = Squonk2Org.objects.all().first()
        if squonk2_org and squonk2_org.uuid != self.__CFG_SQUONK2_ORG_UUID:
            msg: str = f'Configured Squonk2 Organisation ({self.__CFG_SQUONK2_ORG_UUID})'\
                       f' does not match pre-existing record ({squonk2_org.uuid})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # OK, so the ORG UUID has not changed.
        # Is it known to the configured AS?
        if not self._get_squonk2_owner_tokens():
            msg = 'Failed to get AS or DM token for organisation owner'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Get the ORG from the AS API.
        # If it knows the org the response will be successful,
        # and we'll also have the Org's name.
        as_o_rv = AsApi.get_organisation(self.__org_owner_as_token,
                                         org_id=self.__CFG_SQUONK2_ORG_UUID)
        if not as_o_rv.success:
            msg = 'Failed to get AS Organisation'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # The org is known to the AS.
        # Get the AS API version (for reference)
        as_v_rv: AsApiRv = AsApi.get_version()
        if not as_v_rv.success:
            msg = 'Failed to get version from AS'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        as_version: str = as_v_rv.msg['version']

        # If there's no Squonk2Org record, create one,
        # recording the ORG ID and the AS we used to verify it exists.
        if not squonk2_org:
            _LOGGER.info('Creating NEW Squonk2Org record for %s.'
                         ' as-url=%s as-org="%s" as-version=%s',
                         self.__CFG_SQUONK2_ORG_UUID,
                         self.__CFG_SQUONK2_ASAPI_URL,
                         as_o_rv.msg['name'],
                         as_version)
            squonk2_org = Squonk2Org(uuid=self.__CFG_SQUONK2_ORG_UUID,
                                     name=as_o_rv.msg['name'],
                                     as_url=self.__CFG_SQUONK2_ASAPI_URL,
                                     as_version=as_version)
            squonk2_org.save()
        else:
            _LOGGER.debug('Squonk2Org already exists for %s - nothing to do',
                          self.__CFG_SQUONK2_ORG_UUID)

        # Keep the record ID for future use.
        self.__org_record = squonk2_org

        # Organisation is known to AS, and it hasn't changed.
        return SuccessRv

    def _delete_as_product(self, product_uuid: str) -> None:
        """Used in error conditions to remove a previously created Product.
        If this fails there's nothing else we can do so we just return regardless.
        """
        _LOGGER.warning('Deleting AS Product %s', product_uuid)

        as_rv: AsApiRv = AsApi.delete_product(self.__org_owner_as_token,
                                              product_id=product_uuid)
        if not as_rv.success:
            _LOGGER.error('Failed to delete AS Product %s', product_uuid)
            return

        _LOGGER.warning('Deleted AS Product %s', product_uuid)

    def _delete_dm_project(self, project_uuid: str) -> None:
        """Used in error conditions to remove a previously created Project.
        If this fails there's nothing else we can do so we just return regardless.
        """
        _LOGGER.warning('Deleting DM Project %s', project_uuid)

        dm_rv: DmApiRv = DmApi.delete_project(self.__org_owner_dm_token,
                                              project_id=project_uuid)
        if not as_rv.success:
            _LOGGER.error('Failed to delete DM Project %s', project_uuid)
            return

        _LOGGER.warning('Deleted DM Project %s', project_uuid)

    def _create_product_and_project(self,
                                    unit: Squonk2Unit,
                                    user_name: str,
                                    session_title: str,
                                    params: CommonParams) -> Squonk2AgentRv:
        """Called if a Product (and Project) needs to be created. If successful,
        this call returns a dictionary in the response "msg" that contains values
        for the keys "sq2_product_uuid" and "sq2_project_uuid".
        """
        assert unit
        assert user_name
        assert session_title
        assert params

        # Create an AS Product.
        name_truncated, name_full = self._build_product_name(user_name, session_title)
        msg: str = f'Creating AS Product "{name_truncated}"...'
        _LOGGER.info(msg)

        as_rv: AsApiRv = AsApi.create_product(self.__org_owner_as_token,
                                              product_name=name_truncated,
                                              unit_id=unit.uuid,
                                              product_type=_SQ2_PRODUCT_TYPE,
                                              flavour=self.__CFG_SQUONK2_PRODUCT_FLAVOUR)
        if not as_rv.success:
            msg = f'Failed to create AS Product ({as_rv.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        product_uuid: str = as_rv.msg['id']
        msg = f'Created AS Product "{product_uuid}"...'
        _LOGGER.info(msg)

        # Create a DM Project
        name_truncated, name_full = self._build_project_name(1, 2)
        msg = f'Creating DM Project "{name_truncated}"...'
        _LOGGER.info(msg)

        dm_rv: DmApiRv = DmApi.create_project(self.__org_owner_dm_token,
                                              project_name=name_truncated,
                                              as_tier_product_id=product_uuid)
        if not dm_rv.success:
            msg = f'Failed to create DM Project ({dm_rv.msg})'
            _LOGGER.error(msg)
            # First delete the AS Product it should have been attached to
            self._delete_as_product(product_uuid)
            # Then leave...
            return Squonk2AgentRv(success=False, msg=msg)

        project_uuid: str = dm_rv.msg["project_id"]
        msg = f'Created DM Project "{project_uuid}"...'
        _LOGGER.info(msg)

        # Add the user as an Editor to the Project
        msg = f'Adding "{user_name}" to DM Project {project_uuid} as Editor...'
        _LOGGER.info(msg)
        dm_rv = DmApi.add_project_editor(self.__org_owner_dm_token,
                                         project_id=project_uuid,
                                         editor=user_name)
        if not dm_rv.success:
            msg = f'Failed to add "{user_name}" to DM Project ({dm_rv.msg})'
            _LOGGER.error(msg)
            # First delete the DM Project amd the corresponding AS Product...
            self._delete_dm_project(project_uuid)
            self._delete_as_product(product_uuid)
            # Then leave...
            return Squonk2AgentRv(success=False, msg=msg)

        msg = f'Added "{user_name} to DM Project {project_uuid} as Editor'
        _LOGGER.info(msg)

        # If the second call fails - delete the object created in the first

        response_msg: Dict[str, Any] = {"sq2_project_uuid": project_uuid,
                                        "sq2_product_uuid": product_uuid}
        return Squonk2AgentRv(success=True, msg=response_msg)

    def _ensure_unit(self, access_id: int) -> Squonk2AgentRv:
        """Gets or creates a Squonk2 Unit based on a customer's "target access string"
        (TAS). If a Unit is created its name will begin with the text "Fragalysis "
        followed by the configured 'SQUONK2_SLUG' (chosen to be unique between all
        Fragalysis instances that share the same Squonk2 service) and then the
        TAS. In DLS the TAS is essentially the "proposal".

        On success the returned message is used to carry the Squonk2 project UUID.
        """
        assert self.__org_record

        target_access_string = self._get_target_access_string(access_id)
        assert target_access_string

        # Now we check and create a Squonk2Unit...
        unit_name_truncated, unit_name_full = self._build_unit_name(target_access_string)
        sq2_unit: Optional[Squonk2Unit] = Squonk2Unit.objects.filter(name=unit_name_full).first()
        if not sq2_unit:
            _LOGGER.info('No existing Unit for "%s"', target_access_string)
            rv: AsApiRv = AsApi.create_unit(self.__org_owner_as_token,
                                            unit_name=unit_name_truncated,
                                            org_id=self.__org_record.uuid,
                                            billing_day=self.__unit_billing_day)
            if not rv.success:
                msg: str = rv.msg['error']
                _LOGGER.error('Failed to create Unit "%s"', target_access_string)
                return Squonk2AgentRv(success=False, msg=msg)

            unit_uuid: str = rv.msg['id']
            sq2_unit = Squonk2Unit(uuid=unit_uuid,
                                   name=unit_name_full,
                                   organisation_id=self.__org_record.id)
            sq2_unit.save()

            _LOGGER.info('Created NEW Unit %s for "%s"', unit_uuid, target_access_string)
        else:
            _LOGGER.debug('Unit %s already exists for "%s" - nothing to do',
                          sq2_unit.uuid,
                          target_access_string)

        return Squonk2AgentRv(success=True, msg=sq2_unit)

    def _ensure_project(self, params: CommonParams) -> Squonk2AgentRv:
        """Gets or creates a Squonk2 Project, used as the destination of files
        and job executions. Each project requires an AS Product
        (tied to the User and Session) and Unit (tied to the Proposal/Project).

        The proposal is expected to be valid for a given user, this method does not
        check whether the user/proposal combination - it assumes that what's been
        given has been checked.

        On success the returned message is used to carry the Squonk2Project record.

        For testing the target and user IDs are permitted to be 0.
        """
        assert params
        assert isinstance(params, CommonParams)

        # A Squonk2Unit must exist for the Target Access String.
        rv: Squonk2AgentRv = self._ensure_unit(params.access_id)
        if not rv.success:
            return rv
        unit: Squonk2Unit = rv.msg

        user_name: str = self._get_user_name(params.user_id)
        session_title: str = self._get_session_title(params.session_id)
        assert user_name
        assert session_title

        name_truncated, name_full = self._build_product_name(user_name, session_title)
        sq2_project: Optional[Squonk2Project] = Squonk2Project.objects.filter(name=name_full).first()
        if not sq2_project:
            msg = f'No existing Squonk2Project for "{name_full}"'
            _LOGGER.info(msg)
            # Need to call upon Squonk2 to create a 'Product'
            # (and corresponding 'Product').
            rv = self._create_product_and_project(unit, user_name, session_title, params)
            if not rv.success:
                msg = f'Failed creating AS Product or DM Project ({rv.msg})'
                _LOGGER.error(msg)
                return rv

            # Now record these new remote objects in a new
            # Squonk2Project record. As it's worked we're given
            # a dictionary with keys "sq2_project_uuid" and "sq2_product_uuid"
            sq2_project = Squonk2Project(uuid=rv.msg['sq2_project_uuid'],
                                         name=name_full,
                                         product_uuid=rv.msg['sq2_product_uuid'],
                                         unit_id=unit.id)
            sq2_project.save()
            msg = f'Created NEW Squonk2Project for "{name_full}"'
            _LOGGER.info(msg)
        else:
            msg = f'Project {sq2_project.uuid} already exists for "{name_full}" - nothing to do'
            _LOGGER.debug(msg)

        return Squonk2AgentRv(success=True, msg=sq2_project)

    @synchronized
    def get_ui_url(self):
        """Returns the UI URL, if configured.
        """
        return self.__CFG_SQUONK2_UI_URL

    @synchronized
    def configured(self) -> Squonk2AgentRv:
        """Returns True if the module appears to be configured,
        i.e. all the environment variables appear to be set.
        """

        # To prevent repeating the checks, all of which are based on
        # static (environment) variables, if we've been here before
        # just return our previous result.
        if self.__configuration_checked:
            return Squonk2AgentRv(success=self.__configured, msg=None)

        self.__configuration_checked = True
        for name, value in self.__dict__.items():
            # All required configuration has a class '__CFG' prefix
            if name.startswith('_Squonk2Agent__CFG_'):
                if value is None:
                    cfg_name: str = name.split('_Squonk2Agent__CFG_')[1]
                    msg = f'{cfg_name} is not set'
                    _LOGGER.error(msg)
                    return Squonk2AgentRv(success=False, msg=msg)

        # If we get here all the required configuration variables are set

        # Is the slug too long?
        # Limited to 10 characters
        if len(self.__CFG_SQUONK2_SLUG) > _MAX_SLUG_LENGTH:
            msg = f'Slug is longer than {_MAX_SLUG_LENGTH} characters'\
                  f' ({self.__CFG_SQUONK2_SLUG})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Extract hostname and realm from the legacy variable
        # i.e. we need 'example.com' and 'xchem'
        # from 'https://example.com/auth/realms/xchem'
        url: ParseResult = urlparse(self.__CFG_OIDC_KEYCLOAK_REALM)
        self.__keycloak_hostname = url.hostname
        self.__keycloak_realm = os.path.split(url.path)[1]

        # Can we translate the billing day to an integer?
        if not self.__CFG_SQUONK2_UNIT_BILLING_DAY.isdigit():
            msg = 'SQUONK2_UNIT_BILLING_DAY is set'\
                  ' but the value is not a number'\
                  f' ({ self.__CFG_SQUONK2_UNIT_BILLING_DAY})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        self.__unit_billing_day = int(self.__CFG_SQUONK2_UNIT_BILLING_DAY)
        if self.__unit_billing_day < 1:
            msg = 'SQUONK2_UNIT_BILLING_DAY cannot be less than 1'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Product tier flavour must be one of a known value.
        # It's stored in the object's self.__product_flavour as an upper-case value
        if not self.__CFG_SQUONK2_PRODUCT_FLAVOUR in _SUPPORTED_PRODUCT_FLAVOURS:
            msg = f'SQUONK2_PRODUCT_FLAVOUR ({self.__CFG_SQUONK2_PRODUCT_FLAVOUR})' \
                  ' is not supported'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Don't verify Squonk2 SSL certificates?
        if self.__SQUONK2_VERIFY_CERTIFICATES and self.__SQUONK2_VERIFY_CERTIFICATES.lower() == 'no':
            self.__verify_certificates = False
            disable_warnings(InsecureRequestWarning)

        # OK - it all looks good.
        # Mark as 'configured'
        self.__configured = True

        return SuccessRv

    @synchronized
    def ping(self) -> Squonk2AgentRv:
        """Returns True if all the Squonk2 installations
        referred to by the URLs respond.

        We also validate that the organisation supplied is known to the Account Server
        by calling on '_pre_flight_checks()'. If the org is known a Squonk2Org record
        is created, if not the ping fails.
        """
        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        if not self.configured():
            msg = 'Not configured'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Check the UI, DM and AS...

        resp: Optional[Response] = None
        url: str = self.__CFG_SQUONK2_UI_URL
        try:
            resp: Response = requests.head(url, verify=self.__verify_certificates)
        except:
            _LOGGER.error('Exception checking UI at %s', url)
        if resp is None or resp.status_code != 200:
            msg = f'Squonk2 UI is not responding from {url}'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        resp = None
        url = f'{self.__CFG_SQUONK2_DMAPI_URL}/api'
        try:
            resp = requests.head(url, verify=self.__verify_certificates)
        except Exception as ex:
            _LOGGER.error('Exception checking DM at %s', url)
        if resp is None or resp.status_code != 308:
            msg = f'Squonk2 DM is not responding from {url}'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        resp = None
        url = f'{self.__CFG_SQUONK2_ASAPI_URL}/api'
        try:
            resp = requests.head(url, verify=self.__verify_certificates)
        except:
            _LOGGER.error('Exception checking AS at %s', url)
        if resp is None or resp.status_code != 308:
            msg = f'Squonk2 AS is not responding from {url}'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # OK so far.
        # Is the configured organisation known to the AS (and has it changed?)
        status, msg = self._pre_flight_checks()
        if not status:
            msg = f'Failed pre-flight checks ({msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Everything's responding if we get here...
        return SuccessRv

    @synchronized
    def can_run_job(self, params: RunJobParams) -> Squonk2AgentRv:
        """Executes a Job on a Squonk2 installation.
        """
        assert params
        assert isinstance(params, RunJobParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Protect against lack of config or connection/setup issues...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Ensure that the user is allowed to use the given access ID
        user: User  = User.objects.filter(id=params.user_id).first()
        assert user
        access_id: str = params.common.access_id
        assert access_id
        proposal_list: List[str] = self.__ispyb_safe_query_set.get_proposals_for_user(user)
        if not access_id in proposal_list:
            msg = f'The user does not have access to {access_id}'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)
            
        return SuccessRv

    @synchronized
    def can_send(self, params: SendParams) -> Squonk2AgentRv:
        """A blocking method that checks whether a user can send files to Squonk2.
        """
        assert params
        assert isinstance(params, SendParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API**MUST**  call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Ensure that the user is allowed to use the access ID
        user: User  = User.objects.filter(id=params.user_id).first()
        assert user
        access_id: str = params.common.access_id
        assert access_id
        proposal_list: List[str] = self.__ispyb_safe_query_set.get_proposals_for_user(user)
        if not access_id in proposal_list:
            msg = f'You cannot access {access_id}'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)
            
        return SuccessRv

    @synchronized
    def send(self, params: SendParams) -> Squonk2AgentRv:
        """A blocking method that takes care of sending a set of files to
        the configured Squonk2 installation.
        """
        assert params
        assert isinstance(params, SendParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API**MUST**  call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Ensure that the user is allowed to use the access ID
        user = None
        access_id: str = params.common.access_id
        assert access_id
        proposal_list: List[str] = self.__ispyb_safe_query_set.get_proposals_for_user(user)
        if not access_id in proposal_list:
            msg = f'The user does not have access to {access_id}'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)
            
        rv_u: Squonk2AgentRv = self._ensure_project(params.common)
        if not rv_u.success:
            msg = 'Failed to create corresponding Squonk2 Project'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        return SuccessRv

    @synchronized
    def ensure_project(self, params: CommonParams) -> Squonk2AgentRv:
        """A blocking method that takes care of the provisioning of the
        required Squonk2 environment. For Fragalysis this entails the
        creation of a 'Squonk2 Project' (which also requires a 'Unit' and 'Product').

        If successful the Corresponding Squonk2Project record is returned as
        the response 'msg' value.
        """
        assert params
        assert isinstance(params, CommonParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API**MUST**  call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        rv_u: Squonk2AgentRv = self._ensure_project(params)
        if not rv_u.success:
            msg = 'Failed to create corresponding Squonk2 Project'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        return rv_u

# A placeholder for the Agent object
_AGENT_SINGLETON: Optional[Squonk2Agent] = None

def get_squonk2_agent() -> Squonk2Agent:
    """Returns a 'singleton'.
    """
    global _AGENT_SINGLETON  # pylint: disable=global-statement

    if _AGENT_SINGLETON:
        return _AGENT_SINGLETON
    _LOGGER.debug("Creating new Squonk2Agent...")
    _AGENT_SINGLETON = Squonk2Agent()
    _LOGGER.debug("Created")

    return _AGENT_SINGLETON
