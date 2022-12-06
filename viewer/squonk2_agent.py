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
import requests
from requests import Response
from wrapt import synchronized

from viewer.models import Squonk2Project, Squonk2Org, Squonk2Unit

_LOGGER: logging.Logger = logging.getLogger(__name__)

# Response value for the agent methods
Squonk2AgentRv: namedtuple = namedtuple('Squonk2AgentRv', ['success', 'msg'])
SuccessRv: Squonk2AgentRv = Squonk2AgentRv(success=True, msg=None)

# Named tuples are used to pass parameters to the agent methods.
# RunJob, used in run_job()
RunJob: namedtuple = namedtuple("RunJob", ["access_token",
                                           "proposal",
                                           "user_id",
                                           "target_id",
                                           "job_spec",
                                           "callback_url"])

# Send, used in send()
# Access ID is the Fragalysis Project record ID
# Session ID is the Fragalysis SessionProject record ID
# Target ID is the Fragalysis Target record ID
# Snapshot ID is the Fragalysis Snapshot record ID
Send: namedtuple = namedtuple("Send", ["token",
                                       "user_id",
                                       "access_id",
                                       "target_id",
                                       "session_id",
                                       "snapshot_id"])


_SUPPORTED_PRODUCT_FLAVOURS: List[str] = ["BRONZE", "SILVER", "GOLD"]

_MAX_SQ2_NAME_LENGTH: int = 80
_MAX_SLUG_LENGTH: int = 10
_SQ2_NAME_PREFIX: str = "Fragalysis"

# How long are auxiliary parts of Squonk2 names allowed to be?
# This is the maximum name length (assumed to be 80) minus the
# name prefix, minus the slug and minus 2 (for spaces)
_MAX_AUX_STRING_LENGTH: int = _MAX_SQ2_NAME_LENGTH \
                              - len(_SQ2_NAME_PREFIX) \
                              - _MAX_SLUG_LENGTH \
                              - 2

# True if the code's in Test Mode
_TEST_MODE: bool = True


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
            os.environ.get('SQUONK2_SLUG')[:_MAX_SLUG_LENGTH]
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
        self.__FALLBACK_PROPOSAL_ID: Optional[str] =\
            os.environ.get('FALLBACK_PROPOSAL_ID')
        self.__FORCE_FALLBACK_PROPOSAL_ID: Optional[str] =\
            os.environ.get('FORCE_FALLBACK_PROPOSAL_ID')
        self.__SQUONK2_VERIFY_CERTIFICATES: Optional[str] = \
            os.environ.get('SQUONK2_VERIFY_CERTIFICATES')

        # The integer billing day, valid if greater than zero
        self.__unit_billing_day: int = 0
        # The product tier, valid if set
        self.__product_flavour: str = ''
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

        self.__owner_token: str = ''
        self.__keycloak_hostname: str = ''
        self.__keycloak_realm: str = ''

    def _get_org_owner_token(self) -> Optional[str]:
        """Gets an access token for the Squonk2 organisation owner.
        This sets the __keycloak_hostname member and also returns the token.
        """
        assert self.__keycloak_hostname

        _LOGGER.info('__keycloak_hostname="%s" __keycloak_realm="%s" client=%s org_owner=%s',
                     self.__keycloak_hostname,
                     self.__keycloak_realm,
                     self.__CFG_OIDC_AS_CLIENT_ID,
                     self.__CFG_OIDC_AS_CLIENT_ID)

        self.__owner_token = Auth.get_access_token(
            keycloak_url="https://" + self.__keycloak_hostname + "/auth",
            keycloak_realm=self.__keycloak_realm,
            keycloak_client_id=self.__CFG_OIDC_AS_CLIENT_ID,
            username=self.__CFG_SQUONK2_ORG_OWNER,
            password=self.__CFG_SQUONK2_ORG_OWNER_PASSWORD,
        )
        if not self.__owner_token:
            _LOGGER.warning('Failed to get access token for AS Organisation owner')
            return None
        # OK if we get here
        return self.__owner_token

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
            return Squonk2AgentRv(success=False, msg=msg)

        # OK, so the ORG UUID has not changed.
        # Is it known to the configured AS?
        if not self._get_org_owner_token():
            msg = 'Failed to get AS token for organisation owner'
            return Squonk2AgentRv(success=False, msg=msg)

        # Get the ORG from the AS API.
        # If it knows the org the response will be successful,
        # and we'll also have the Org's name.
        _LOGGER.info('Checking organisation (%s)', self.__CFG_SQUONK2_ORG_UUID)
        as_o_rv = AsApi.get_organisation(self.__owner_token,
                                         org_id=self.__CFG_SQUONK2_ORG_UUID)
        if not as_o_rv.success:
            msg = 'Failed checking AS Organisation'
            return Squonk2AgentRv(success=False, msg=msg)

        # The org is known to the AS.
        # Get the AS API version (for reference)
        as_v_rv: AsApiRv = AsApi.get_version()
        if not as_v_rv.success:
            msg = 'Failed to get version from AS'
            return Squonk2AgentRv(success=False, msg=msg)
        as_version: str = as_v_rv.msg['version']

        # If there's no Squonk2Org record, create one,
        # recording the ORG ID and the AS we used to verify it exists.
        if not squonk2_org:
            squonk2_org = Squonk2Org(uuid=self.__CFG_SQUONK2_ORG_UUID,
                                     name=as_o_rv.msg['name'],
                                     as_url=self.__CFG_SQUONK2_ASAPI_URL,
                                     as_version=as_version)
            squonk2_org.save()

        # Keep the record ID for future use.
        self.__org_record = squonk2_org

        # Organisation is known to AS, and it hasn't changed.
        return SuccessRv

    def _ensure_unit(self, access_id: int) -> Squonk2AgentRv:
        """Gets or creates a Squonk2 Unit based on a customer's "target access string"
        (TAS). If a Unit is created its name will begin with the text "Fragalysis "
        followed by the configured 'SQUONK2_SLUG' (chosen to be unique between all
        Fragalysis instances that share the same Squonk2 service) and then the
        TAS. In DLS the TAS is essentially the "proposal".

        On success the returned message is used to carry the Squonk2 project UUID.
        """
        assert self.__org_record
        if not _TEST_MODE:
            assert access_id

        # Get the target access string
        target_access_string: str = ''
        if _TEST_MODE:
            _LOGGER.warning('Using FALLBACK_PROPOSAL_ID')
            target_access_string = self.__FALLBACK_PROPOSAL_ID

        unit: Optional[Squonk2Unit] = Squonk2Unit.objects.filter(proposal=target_access_string).first()
        if not unit:
            _LOGGER.info('No existing Unit for this TAS (%s)', target_access_string)
            unit_name: styr = 'Fragalysis' \
                              f' {self.__CFG_SQUONK2_SLUG}' \
                              f' {target_access_string[:_MAX_AUX_STRING_LENGTH]}'
            rv: AsApiRv = AsApi.create_unit(unit_name=unit_name,
                                            org_id=self.__org_record.uuid,
                                            billing_day=self.__unit_billing_day)
            if not rv.success:
                msg: str = rv.msg['error']
                _LOGGER.error('Failed to create Unit for TAS (%s)', target_access_string)
                return Squonk2AgentRv(success=False, msg=msg)

            unit_uuid: str = rv.msg['id']
            unit: Squonk2Unit = Squonk2Unit(uuid=unit_uuid,
                                            name=unit_name,
                                            proposal=target_access_string,
                                            organisation=self.__org_record.id)
            unit.save()
        else:
            _LOGGER.debug('Unit %s already exists for this TAS (%s)',
                          unit.uuid,
                          target_access_string)

        return Squonk2AgentRv(success=True, msg=unit.uuid)

    def _ensure_project(self,
                        target_id: int,
                        user_id: int) -> Squonk2AgentRv:
        """Gets or creates a Squonk2 Project, used as the destination of files
        and job executions. Each project requires an AS Product
        (tied to the User and Target) and Unit (tied to the proposal).

        The proposal is expected to be valid for a given user, this method does not
        check whether the user/proposal combination - it assumes that what's been
        given has been checked.

        On success the returned message is used to carry the Squonk2 project UUID.

        For testing the target and user IDs are permitted to be 0.
        """
        if not _TEST_MODE:
            assert target_id > 0
            assert user_id > 0

        # A Squonk2Unit must exist for the Target Access String, and there must be
        # a Squonk2Project record for the user/target combination.
        # If not they are created.

        rv: Squonk2AgentRv = self._ensure_unit(target_id)
        if not rv.success:
            return rv
        unit_uuid: str = rv.msg
        print(f'Got or Created Unit {unit_uuid}')

        # A Squonk2Project record must exist for the unit/user/target combination.
        # If not it is created.
        project: Optional[Squonk2Project] =\
            Squonk2Project.objects.filter(user__id=user_id,
                                          target__id=target_id,
                                          unit__uuid=unit_uuid)\
                .first()
        if not project:
            # Need to call upon Squonk2 to create a 'Product' and a 'Project'.
            # The Product is created by the organisation 'owner'
            # the 'Project' is created on behalf of the 'user'.
            #
            # We create a corresponding Squonk2Project when both have been successful.
            rv = self._create_project()
            if not rv.success:
                return rv
            project = Squonk2Project(uuid=project_uuid,
                                     name=project_name,
                                     product_uuid=product_uuid,
                                     unit=unit_id,
                                     user=user_id,
                                     target=target_id)
            project.save()

        return Squonk2AgentRv(success=True, msg=project.uuid)

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
                    msg: str = f'{cfg_name} is not set'
                    _LOGGER.warning(msg)
                    return Squonk2AgentRv(success=False, msg=msg)

        # If we get here all the required configuration variables are set

        # Is the slug too long?
        # Limited to 10 characters
        if len(self.__CFG_SQUONK2_SLUG) > _MAX_SLUG_LENGTH:
            msg: str = f'Slug is longer than {_MAX_SLUG_LENGTH} characters'\
                       f' ({self.__CFG_SQUONK2_SLUG})'
            _LOGGER.warning(msg)
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

        # Product tier to upper-case
        if not self.__CFG_SQUONK2_PRODUCT_FLAVOUR in _SUPPORTED_PRODUCT_FLAVOURS:
            msg = f'SQUONK2_PRODUCT_FLAVOUR ({self.__CFG_SQUONK2_PRODUCT_FLAVOUR})' \
                  ' is not supported'
            _LOGGER.error(msg)
            return Squonk2AgentRv(success=False, msg=msg)
        self.__product_flavour = self.__CFG_SQUONK2_PRODUCT_FLAVOUR.upper()

        # Don't verify Squonk2 SSL certificates?
        if self.__SQUONK2_VERIFY_CERTIFICATES and self.__SQUONK2_VERIFY_CERTIFICATES.lower() == 'no':
            self.__verify_certificates = False
            disable_warnings(InsecureRequestWarning)

        # OK - it all looks good.
        # Mark as 'configured'
        self.__configured = True

        return SuccessRv

    def ping(self) -> Squonk2AgentRv:
        """Returns True if all the Squonk2 installations
        referred to by the URLs respond.

        We also validate that the organisation supplied is known to the Account Server
        by calling on '_pre_flight_checks()'. If the org is known a Squonk2Org record
        is created, if not the ping fails.
        """
        if not self.configured():
            msg: str = 'Not configured'
            _LOGGER.debug(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # Check the UI, DM and AS...

        resp: Optional[Response] = None
        url: str = self.__CFG_SQUONK2_UI_URL
        try:
            resp: Response = requests.head(url, verify=self.__verify_certificates)
        except:
            _LOGGER.warning('Exception checking UI at %s', url)
        if resp is None or resp.status_code != 200:
            msg = f'UI is not responding from {url}'
            _LOGGER.debug(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        resp = None
        url = f'{self.__CFG_SQUONK2_DMAPI_URL}/api'
        try:
            resp = requests.head(url, verify=self.__verify_certificates)
        except Exception as ex:
            _LOGGER.warning('Exception checking DM at %s', url)
        if resp is None or resp.status_code != 308:
            msg = f'Data Manager is not responding from {url}'
            _LOGGER.debug(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        resp = None
        url = f'{self.__CFG_SQUONK2_ASAPI_URL}/api'
        try:
            resp = requests.head(url, verify=self.__verify_certificates)
        except:
            _LOGGER.warning('Exception checking AS at %s', url)
        if resp is None or resp.status_code != 308:
            msg = f'Account Manager is not responding from {url}'
            _LOGGER.debug(msg)
            return Squonk2AgentRv(success=False, msg=msg)

        # OK so far.
        # Is the configured organisation known to the AS (and has it changed?)
        status, msg = self._pre_flight_checks()
        if not status:
            msg = f'Failed pre-flight checks ({msg})'
            return Squonk2AgentRv(success=False, msg=msg)

        # Everything's responding if we get here...
        return SuccessRv

    def run_job(self, params: RunJob) -> Squonk2AgentRv:
        """Executes a Job on a Squonk2 installation.
        """
        assert params
        assert isinstance(params, RunJob)

        # Protect against lack of config or connection/setup issues...
        if not self.ping():
            msg: str = 'Squonk2 ping failed.'\
                       ' Are we configured properly and is Squonk alive?'
            return Squonk2AgentRv(success=False, msg=msg)

        return SuccessRv

    def send(self, params: Send) -> Squonk2AgentRv:
        """A blocking method that takes care of sending a set of files to
        the configured Squonk2 installation.
        """
        assert params
        assert isinstance(params, Send)

        # Protect against lack of config or connection/setup issues...
        if not self.ping():
            msg: str = 'Squonk2 ping failed.'\
                       ' Are we configured properly and is Squonk alive?'
            return Squonk2AgentRv(success=False, msg=msg)

        rv_u: Squonk2AgentRv = self._ensure_project(params.user_id,
                                                    params.access_id,
                                                    params.target_id)
        if not rv_u.succes:
            msg: str = 'Failed to create corresponding Squonk2 Project'
            return Squonk2AgentRv(success=False, msg=msg)

        return SuccessRv


# The global (singleton).
# This acts as out sole singleton,
# created and returned from `get_squonk2_agent()'
_SQUONK2_AGENT: Optional[Squonk2Agent] = None


def get_squonk2_agent() -> Squonk2Agent:
    """Returns a 'singleton'.
    """
    global _SQUONK2_AGENT  # pylint: disable=global-statement

    if _SQUONK2_AGENT:
        return _SQUONK2_AGENT
    _LOGGER.debug("Creating new Squonk2Agent...")
    _SQUONK2_AGENT = Squonk2Agent()
    _LOGGER.debug("Created")

    return _SQUONK2_AGENT
