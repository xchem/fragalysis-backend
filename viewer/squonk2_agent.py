"""An agent for the Squonk2 (Data Manger and Account Server) API.
This module 'simplifies' the use of the Squonk2 Python client package.
"""

# Refer to the accompanying low-level-design document: -
# https://docs.google.com/document/d/1lFpN29dK1luz80lwSGi0Rnj1Rqula2_CTRuyWDUBu14

from collections import namedtuple
import logging
import os
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import ParseResult, urlparse
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings

from squonk2.auth import Auth
from squonk2.as_api import AsApi, AsApiRv
from squonk2.dm_api import DmApi, DmApiRv
import requests
from requests import Response
from wrapt import synchronized

from api.security import ISpyBSafeQuerySet
from viewer.models import User, Project, Target
from viewer.models import Squonk2Project, Squonk2Org, Squonk2Unit

_LOGGER: logging.Logger = logging.getLogger(__name__)

# Response value for the agent methods.
# It's a boolean and an optional string (used for errors or response content)
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

# Parameters for the grant_access() method.
AccessParams: namedtuple = namedtuple("Access", ["username",
                                                 "project_uuid"])

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
        self.__DUMMY_TARGET_TITLE: Optional[str] =\
            os.environ.get('DUMMY_TARGET_TITLE')
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
        # Ignore cert errors? (no)
        self.__verify_certificates: bool = True

        # The record ID of the Squonk2Org for this deployment.
        # Set on successful 'pre-flight-check'
        self.__org_record: Optional[Squonk2Org] = None

        self.__org_owner_as_token: str = ''
        self.__org_owner_dm_token: str = ''
        self.__keycloak_hostname: str = ''
        self.__keycloak_realm: str = ''

        # The Safe QuerySet from the security module.
        # Used when we are given a tas (target access string).
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

    def _get_target_title(self, target_id: int) -> str:
        # Gets the Target title (if it looks sensible)
        # or a fixed value (used for testing)
        if target_id == 0:
            assert _TEST_MODE
            _LOGGER.warning('Caution - in TEST mode, using __DUMMY_TARGET_TITLE (%s)',
                            self.__DUMMY_TARGET_TITLE)

        target_title: str = self.__DUMMY_TARGET_TITLE
        if target_id:
            target: Target = Target.objects.filter(id=target_id).first()
            assert target
            target_title = target.title
        assert target_title
        return target_title

    def _build_unit_name(self, target_access_string: str) -> Tuple[str, str]:
        assert target_access_string
        # AS Units are named using the Target Access String (TAS)
        # which, in Diamond will be a "visit" string like "lb00000-1"
        name: str = f'{_SQ2_NAME_PREFIX} {self.__CFG_SQUONK2_SLUG} /{target_access_string}/'
        return name[:_SQ2_MAX_NAME_LENGTH], name

    def _build_product_name(self, username: str, target_title: str) -> Tuple[str, str]:
        """Builds a Product name, returning the truncated and un-truncated form"""
        assert username
        assert target_title
        # AS Products are named using the user and the session
        # (there's a 1:1 mapping to DM Projects)

        # The Product name characters are not restricted
        identifier: str = f'{username}::{target_title}'
        name: str = f'{_SQ2_NAME_PREFIX} {self.__CFG_SQUONK2_SLUG} {identifier}'
        return name[:_SQ2_MAX_NAME_LENGTH], name

    def _get_squonk2_owner_tokens(self) -> Optional[Tuple[str, str]]:
        """Gets access tokens for the Squonk2 organisation owner.
        This sets the __keycloak_hostname member and also returns the tokens,
        getting one for the AS and one for the DM.
        """
        assert self.__keycloak_hostname

        _LOGGER.debug('__keycloak_hostname="%s" __keycloak_realm="%s"'
                      ' dm-client=%s as-client=%s org=%s org_owner=%s',
                     self.__keycloak_hostname,
                     self.__keycloak_realm,
                     self.__CFG_OIDC_DM_CLIENT_ID,
                     self.__CFG_OIDC_AS_CLIENT_ID,
                     self.__CFG_SQUONK2_ORG_UUID,
                     self.__CFG_SQUONK2_ORG_OWNER)

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

        assert self.__configuration_checked
        assert self.__configured

        if self.__org_record and self.__org_record.uuid != self.__CFG_SQUONK2_ORG_UUID:
            msg: str = f'Configured Squonk2 Organisation ({self.__CFG_SQUONK2_ORG_UUID})'\
                       f' does not match pre-existing record ({self.__org_record.uuid})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # OK, so the ORG exists and its UUID has not changed.
        # Is it known to the configured AS?
        if not self._get_squonk2_owner_tokens():
            msg = 'Failed to get AS or DM token for organisation owner'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        _LOGGER.debug('Got Squonk2 API Access Tokens')

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
        _LOGGER.debug('Happy with Squonk2 Account Server (as_version=%s)', as_version)

        # Everything seems to be OK but we might not have an organisation in this
        # call (it may be the first call of the instance lifetime).
        # So, if there's no Squonk2Org record, create one,
        # recording the ORG ID and the AS and version we used.
        if not self.__org_record:
            assert self.__CFG_SQUONK2_ASAPI_URL
            _LOGGER.info('Creating NEW Squonk2Org record for %s.'
                         ' as-url=%s as-org="%s" as-version=%s',
                         self.__CFG_SQUONK2_ORG_UUID,
                         self.__CFG_SQUONK2_ASAPI_URL,
                         as_o_rv.msg['name'],
                         as_version)
            self.__org_record = Squonk2Org(uuid=self.__CFG_SQUONK2_ORG_UUID,
                                           name=as_o_rv.msg['name'],
                                           as_url=self.__CFG_SQUONK2_ASAPI_URL,
                                           as_version=as_version)
            self.__org_record.save()
            _LOGGER.info('Created Squonk2Org record for %s',
                         self.__CFG_SQUONK2_ORG_UUID)
        else:
            _LOGGER.debug('Squonk2Org for %s "%s" already exists - nothing to do',
                          self.__org_record.uuid,
                          self.__org_record.name)

        # Organisation is known to AS, and it hasn't changed.
        _LOGGER.debug('Successful pre-flight checks')
        return SuccessRv

    def _get_or_create_unit(self,
                            unit_name_truncated: str,
                            unit_name_full: str) -> Squonk2AgentRv:
        """Gets an exiting Unit or creates a new one
        returning its UUID as the msg content.
        """
        # Get existing Units for our Organisation
        org_uuid: str = self.__org_record.uuid
        as_rv: AsApiRv = AsApi.get_units(self.__org_owner_as_token, org_id=org_uuid)
        if not as_rv.success:
            msg: str = as_rv.msg['error']
            _LOGGER.error('Failed to get Units for Organisation "%s"', org_uuid)
            return Squonk2AgentRv(success=False, msg=msg)
        # Iterate through them, looking for ours...
        _LOGGER.info('Got %s Units for Organisation "%s"', len(as_rv.msg['units']), org_uuid)
        for unit in as_rv.msg['units']:
            if unit['name'] == unit_name_truncated:
                unit_uuid: str = unit['id']
                _LOGGER.info('...and one was ours (%s)', unit_uuid)
                return Squonk2AgentRv(success=True, msg=unit_uuid)

        # Not found, create it...
        _LOGGER.info('No exiting Unit, creating NEW Unit "%s" (for "%s")',
                      unit_name_full,
                      unit_name_truncated)
        as_rv = AsApi.create_unit(self.__org_owner_as_token,
                                  unit_name=unit_name_truncated,
                                  org_id=org_uuid,
                                  billing_day=self.__unit_billing_day)
        if not as_rv.success:
            msg = as_rv.msg['error']
            _LOGGER.error('Failed to create Unit "%s" (for "%s")',
                          unit_name_full, unit_name_truncated)
            return Squonk2AgentRv(success=False, msg=msg)

        unit_uuid = as_rv.msg['id']
        _LOGGER.info('Created NEW Unit "%s"', unit_uuid)
        return Squonk2AgentRv(success=True, msg=unit_uuid)

    def _get_or_create_product(self, name_truncated: str, unit_uuid: str) -> Squonk2AgentRv:

        # Get existing Products for the Unit
        as_rv: AsApiRv = AsApi.get_products_for_unit(self.__org_owner_as_token,
                                                     unit_id=unit_uuid)
        if not as_rv.success:
            msg: str = as_rv.msg['error']
            _LOGGER.error('Failed to get Products for Unit "%s"', unit_uuid)
            return Squonk2AgentRv(success=False, msg=msg)

        # Iterate through them, looking for ours...
        for product in as_rv.msg['products']:
            if product['product']['name'] == name_truncated:
                product_uuid: str = product['product']['id']
                _LOGGER.info('Found pre-existing AS Product "%s"', product_uuid)
                return Squonk2AgentRv(success=True, msg=product_uuid)

        # No existing Product with this name...        
        _LOGGER.info('No existing Product, creating NEW AS Product "%s" (for "%s")',
                      name_truncated,
                      unit_uuid)
        as_rv = AsApi.create_product(self.__org_owner_as_token,
                                     product_name=name_truncated,
                                     unit_id=unit_uuid,
                                     product_type=_SQ2_PRODUCT_TYPE,
                                     flavour=self.__CFG_SQUONK2_PRODUCT_FLAVOUR)
        if not as_rv.success:
            msg = f'Failed to create AS Product ({as_rv.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        product_uuid = as_rv.msg['id']
        msg = f'Created NEW AS Product {product_uuid}...'
        _LOGGER.info(msg)
        return Squonk2AgentRv(success=True, msg=product_uuid)

    def _get_or_create_project(self, name_truncated: str, product_uuid: str) -> Squonk2AgentRv:
        """Gets existing DM Projects (that belong to the given Product)
        to see if ours exists, if not a new one is created.
        """
        # TODO
        dm_rv: DmApiRv = DmApi.get_available_projects(self.__org_owner_as_token)
        if not dm_rv.success:
            msg: str = dm_rv.msg['error']
            _LOGGER.error('Failed to get Projects')
            return Squonk2AgentRv(success=False, msg=msg)

        # Iterate through them, looking for ours...
        for project in dm_rv.msg['projects']:
            if project['name'] == name_truncated and 'product_id' in project and project['product_id'] == product_uuid:
                project_uuid: str = project['project_id']
                _LOGGER.info('Found pre-existing DM Project "%s"', project_uuid)
                return Squonk2AgentRv(success=True, msg=project_uuid)

        # No existing Project
        _LOGGER.info('No existing Project, creating NEW DM Project "%s" (for "%s")',
                      name_truncated,
                      product_uuid)
        dm_rv = DmApi.create_project(self.__org_owner_dm_token,
                                     project_name=name_truncated,
                                     private=True,
                                     as_tier_product_id=product_uuid)
        if not dm_rv.success:
            msg = f'Failed to create DM Project ({dm_rv.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        
        project_uuid: str = dm_rv.msg['project_id']
        msg = f'Created NEW DM Project {project_uuid}...'
        return Squonk2AgentRv(success=True, msg=project_uuid)

    def _delete_as_product(self, product_uuid: str) -> None:
        """Used in error conditions to remove a previously created Product.
        If this fails there's nothing else we can do so we just return regardless.
        """
        _LOGGER.warning('Deleting AS Product %s...', product_uuid)

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
        _LOGGER.warning('Deleting DM Project %s...', project_uuid)

        dm_rv: DmApiRv = DmApi.delete_project(self.__org_owner_dm_token,
                                              project_id=project_uuid)
        if not dm_rv.success:
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
        name_truncated, _ = self._build_product_name(user_name, session_title)
        msg: str = f'Creating AS Product "{name_truncated}" (unit={unit.uuid})...'
        _LOGGER.info(msg)

        sq2_rv: Squonk2AgentRv = self._get_or_create_product(name_truncated, unit.uuid)
        if not sq2_rv.success:
            msg = f'Failed to create AS Product ({sq2_rv.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        product_uuid: str = sq2_rv.msg
        msg = f'Got or created AS Product {product_uuid}'
        _LOGGER.info(msg)

        # Create a DM Project (using the same name we used for the AS Product)
        msg = f'Continuing by creating DM Project "{name_truncated}"...'
        _LOGGER.info(msg)

        sq2_rv: Squonk2AgentRv = self._get_or_create_project(name_truncated, product_uuid)
        if not sq2_rv.success:
            msg = f'Failed to create DM Project ({sq2_rv.msg})'
            _LOGGER.error(msg)
            # First delete the AS Product it should have been attached to
            self._delete_as_product(product_uuid)
            # Then leave...
            return Squonk2AgentRv(success=False, msg=msg)
        project_uuid: str = sq2_rv.msg
        msg = f'Got or created DM Project {project_uuid}...'
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
            _LOGGER.warning('Rolling back DM Project and AS Product creation...')
            # First delete the DM Project amd the corresponding AS Product...
            self._delete_dm_project(project_uuid)
            self._delete_as_product(product_uuid)
            # Then leave...
            return Squonk2AgentRv(success=False, msg=msg)

        msg = f'"{user_name}" is now an Editor of DM Project {project_uuid}'
        _LOGGER.info(msg)

        # If the second call fails - delete the object created in the first

        response_msg: Dict[str, Any] = {"sq2_project_uuid": project_uuid,
                                        "sq2_product_uuid": product_uuid}
        return Squonk2AgentRv(success=True, msg=response_msg)

    def _ensure_unit(self, target_access_string: int) -> Squonk2AgentRv:
        """Gets or creates a Squonk2 Unit based on a customer's "target access string"
        (TAS). If a Unit is created its name will begin with the text "Fragalysis "
        followed by the configured 'SQUONK2_SLUG' (chosen to be unique between all
        Fragalysis instances that share the same Squonk2 service) and then the
        TAS. In DLS the TAS is essentially the "proposal".

        On success the returned message is used to carry the Squonk2 project UUID.
        """
        if not self.__org_record:
            msg: str = 'The Squonk2Org record does not match' \
                       ' the configured SQUONK2_ORG_UUID.' \
                       ' You cannot change the SQUONK2_ORG_UUID once it has been used'
            return Squonk2AgentRv(success=False, msg=msg)

        # Now we check and create a Squonk2Unit...
        unit_name_truncated, unit_name_full = self._build_unit_name(target_access_string)
        sq2_unit: Optional[Squonk2Unit] = Squonk2Unit.objects.filter(name=unit_name_full).first()
        if not sq2_unit:
            _LOGGER.info('No existing Squonk2Unit for "%s"', target_access_string)
            # Get the list of Units from Squonk.
            sq2a_rv: Squonk2AgentRv = self._get_or_create_unit(unit_name_truncated, unit_name_full)
            if not sq2a_rv.success:
                _LOGGER.error('Failed to create Unit "%s" (%s)', target_access_string, sq2a_rv.msg)
                return Squonk2AgentRv(success=False, msg=sq2a_rv.msg)

            unit_uuid: str = sq2a_rv.msg
            sq2_unit = Squonk2Unit(uuid=unit_uuid,
                                   name=unit_name_full,
                                   organisation_id=self.__org_record.id)
            sq2_unit.save()
            _LOGGER.info('Created Squonk2Unit %s "%s" (for "%s")',
                         unit_uuid,
                         unit_name_full,
                         target_access_string)
        else:
            _LOGGER.debug('Squonk2Unit %s "%s" already exists (for "%s") - nothing to do',
                          sq2_unit.uuid,
                          unit_name_full,
                          target_access_string)

        return Squonk2AgentRv(success=True, msg=sq2_unit)

    def _ensure_project(self, c_params: CommonParams) -> Squonk2AgentRv:
        """Gets or creates a Squonk2 Project, used as the destination of files
        and job executions. Each Project requires an AS Product
        (tied to the User and Session) and Unit (tied to the Proposal/Project).

        The proposal is expected to be valid for a given user, this method does not
        check whether the user/proposal combination - it assumes that what's been
        given has been checked.

        On success the returned message is used to carry the Squonk2Project record.

        For testing the target and user IDs are permitted to be 0.
        """
        assert c_params
        assert isinstance(c_params, CommonParams)
        
        target_access_string = self._get_target_access_string(c_params.access_id)
        assert target_access_string

        # A Squonk2Unit must exist for the Target Access String.
        rv: Squonk2AgentRv = self._ensure_unit(target_access_string)
        if not rv.success:
            return rv
        unit: Squonk2Unit = rv.msg

        user_name: str = self._get_user_name(c_params.user_id)
        target_title: str = self._get_target_title(c_params.target_id)
        assert user_name
        assert target_title

        _, name_full = self._build_product_name(user_name, target_title)
        sq2_project: Optional[Squonk2Project] = Squonk2Project.objects.filter(name=name_full).first()
        if not sq2_project:
            msg = f'No existing Squonk2Project for "{name_full}"'
            _LOGGER.info(msg)
            # Need to call upon Squonk2 to create a 'Product'
            # (and corresponding 'Product').
            rv = self._create_product_and_project(unit, user_name, target_title, c_params)
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
            msg = f'Created NEW Squonk2Project for {sq2_project.uuid} "{name_full}"'
            _LOGGER.info(msg)
        else:
            msg = f'Squonk2Project for {sq2_project.uuid} "{name_full}" already exists - nothing to do'
            _LOGGER.debug(msg)

        return Squonk2AgentRv(success=True, msg=sq2_project)

    def _verify_access(self, c_params: CommonParams) -> Squonk2AgentRv:
        """Checks the user has access to the project.
        """
        access_id: int = c_params.access_id
        assert access_id
        project: Optional[Project] = Project.objects.filter(id=access_id).first()
        if not project:
            msg = f'Access ID (Project) {access_id} does not exist'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Ensure that the user is allowed to use the given access ID.
        # Even on public projects the user must be part of the project
        # to use Squonk.
        user: User  = User.objects.filter(id=c_params.user_id).first()
        assert user
        target_access_string = self._get_target_access_string(access_id)
        assert target_access_string
        proposal_list: List[str] = self.__ispyb_safe_query_set.get_proposals_for_user(user)
        if not target_access_string in proposal_list:
            msg = f'The user ({user.username}) cannot access "{target_access_string}"' \
                  f' (access_id={access_id}). Only {proposal_list})'
            _LOGGER.warning(msg)
            return Squonk2AgentRv(success=False, msg=msg)
            
        return SuccessRv
        
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
            _LOGGER.debug('Configuration already checked (configured=%s)',
                          self.__configured)
            return Squonk2AgentRv(success=self.__configured, msg=None)

        self.__configuration_checked = True
        self.__configured = False

        _LOGGER.debug('Checking configuration (unchecked)...')

        for name, value in self.__dict__.items():
            # All required configuration has a class '__CFG' prefix
            if name.startswith('_Squonk2Agent__CFG_'):
                _LOGGER.debug('CFG variable "%s" has value "%s"', name, value)
                if value is None or not value:
                    cfg_name: str = name.split('_Squonk2Agent__CFG_')[1]
                    msg = f'{cfg_name} is not set'
                    _LOGGER.error(msg)
                    return Squonk2AgentRv(success=False, msg=msg)

        # If we get here all the required configuration variables are set
        _LOGGER.debug('Required configuration variables are set')

        # Is the slug too long?
        # Limited to 10 characters
        if len(self.__CFG_SQUONK2_SLUG) > _MAX_SLUG_LENGTH:
            msg = f'Slug is longer than {_MAX_SLUG_LENGTH} characters'\
                  f' ({self.__CFG_SQUONK2_SLUG})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        _LOGGER.debug('Acceptable slug (%s)', self.__CFG_SQUONK2_SLUG)

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
        _LOGGER.debug('Acceptable billing day (%s)', self.__unit_billing_day)

        # Product tier flavour must be one of a known value.
        # It's stored in the object's self.__product_flavour as an upper-case value
        if not self.__CFG_SQUONK2_PRODUCT_FLAVOUR in _SUPPORTED_PRODUCT_FLAVOURS:
            msg = f'SQUONK2_PRODUCT_FLAVOUR ({self.__CFG_SQUONK2_PRODUCT_FLAVOUR})' \
                  ' is not supported'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        _LOGGER.debug('Acceptable flavour (%s)', self.__CFG_SQUONK2_PRODUCT_FLAVOUR)

        # Don't verify Squonk2 SSL certificates?
        if self.__SQUONK2_VERIFY_CERTIFICATES and self.__SQUONK2_VERIFY_CERTIFICATES.lower() == 'no':
            self.__verify_certificates = False
            disable_warnings(InsecureRequestWarning)
            _LOGGER.debug('Disabled certificate verification')

        # OK - it all looks good.
        # This is the only place where we set '__configured'
        self.__configured = True
        _LOGGER.debug('Configuration checked (configured=%s)', self.__configured)

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
        except Exception:  # pylint: disable=broad-except
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
    def can_run_job(self, rj_params: RunJobParams) -> Squonk2AgentRv:
        """Executes a Job on a Squonk2 installation.
        """
        assert rj_params
        assert isinstance(rj_params, RunJobParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Protect against lack of config or connection/setup issues...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        return self._verify_access(c_params=rj_params.common)

    @synchronized
    def can_send(self, s_params: SendParams) -> Squonk2AgentRv:
        """A blocking method that checks whether a user can send files to Squonk2.
        """
        assert s_params
        assert isinstance(s_params, SendParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API **MUST** call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        return self._verify_access(c_params=s_params.common)

    @synchronized
    def send(self, s_params: SendParams) -> Squonk2AgentRv:
        """A blocking method that takes care of sending a set of files to
        the configured Squonk2 installation.
        """
        assert s_params
        assert isinstance(s_params, SendParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API **MUST** call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        rv_access: Squonk2AgentRv = self._verify_access(s_params.common)
        if not rv_access.success:
            return rv_access
            
        rv_u: Squonk2AgentRv = self._ensure_project(s_params.common)
        if not rv_u.success:
            msg = f'Failed to create corresponding Squonk2 Project (msg={rv_u.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        return SuccessRv

    @synchronized
    def ensure_project(self, c_params: CommonParams) -> Squonk2AgentRv:
        """A blocking method that takes care of the provisioning of the
        required Squonk2 environment. For Fragalysis this entails the
        creation of a 'Squonk2 Project' (which also requires a 'Unit' and 'Product').

        If successful the Corresponding Squonk2Project record is returned as
        the response 'msg' value.
        """
        assert c_params
        assert isinstance(c_params, CommonParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API **MUST** call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        rv_access: Squonk2AgentRv = self._verify_access(c_params)
        if not rv_access.success:
            return rv_access
                        
        rv_u: Squonk2AgentRv = self._ensure_project(c_params)
        if not rv_u.success:
            msg = f'Failed to create corresponding Squonk2 Project (msg={rv_u.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        return rv_u

    @synchronized
    def grant_access(self, a_params: AccessParams) -> Squonk2AgentRv:
        """A blocking method that takes care of sending a set of files to
        the configured Squonk2 installation.
        """
        assert a_params
        assert isinstance(a_params, AccessParams)

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API **MUST**  call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        dm_rv: DmApiRv = DmApi.add_project_observer(self.__org_owner_dm_token,
                                                    project_id=a_params.project_uuid,
                                                    observer=a_params.username)
        if not dm_rv.success:
            msg = f'Failed to add DM Project Observer ({dm_rv.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        
        return SuccessRv

    @synchronized
    def get_instance_execution_status(self, callback_context: str) -> Squonk2AgentRv:
        """A blocking method that attempt to get the execution status (success/failure)
        of an instance (a Job) based on the given callback context. The status (string)
        is returned as the Squonk2AgentRv.msg value.
        """
        assert callback_context

        if _TEST_MODE:
            msg: str = 'Squonk2Agent is in TEST mode'
            _LOGGER.warning(msg)

        # Every public API **MUST**  call ping().
        # This ensures Squonk2 is available and gets suitable API tokens...
        if not self.ping():
            msg = 'Squonk2 ping failed.'\
                  ' Are we configured properly and is Squonk2 alive?'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # To do this we actually get the DM tasks associated with the
        # callback context that Fragalysis provided.
        # For the filters we provide we should only get one Task.
        dm_rv: DmApiRv = DmApi.get_tasks(self.__org_owner_dm_token,
                                         exclude_removal=True,
                                         exclude_purpose='FILE.DATASET.PROJECT',
                                         instance_callback_context=callback_context)
        if not dm_rv.success:
            msg = f'Failed to get DM Tasks ({dm_rv.msg})'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        
        # Find the task that relates to our instance.
        # We return 'LOST' if a task cannot be found,
        # otherwise it's one of None, SUCCESS or FAILURE
        i_status: Optional[str] = 'LOST'
        num_tasks: int = len(dm_rv.msg['tasks'])
        if num_tasks == 1:
            i_task = dm_rv.msg['tasks'][0]
            if i_task['done']:
                i_status = 'FAILURE' if i_task['exit_code'] != 0 else 'SUCCESS'
            else:
                i_status = None
        else:
            msg = f'More than one Task found ({num_tasks}) for callback context "{callback_context}"'
            _LOGGER.warning(msg)

        if i_status and i_status == 'LOST':
            msg = f'No Task found for callback context "{callback_context}", assume "LOST"'
            _LOGGER.warning(msg)
                
        return Squonk2AgentRv(success=True, msg=i_status)

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
