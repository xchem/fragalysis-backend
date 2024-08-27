"""Django settings for the fragalysis 'backend'"""

# This standard Django module is used to provide the dynamic configuration of
# the backend logic. As well as providing vital django-related configuration
# it is also the source of the numerous fragalysis-specific environment variables
# that control the stack's configuration (behaviour).
#
# Not all settings are configured by environment variable. Some are hard-coded
# and you'll need to edit their values here. For example `ALLOWED_HOSTS`
# is a static variable that is not intended to be changed at run-time.
#
# Those that are configurable at run-time should be obvious
# (i.e. they'll use "os.environ.get()" to obtain their value) alternative run-time value.
#
# You will find the django-related configuration at the top of the file
# (under DJANGO SETTINGS) and the fragalysis-specific configuration at the bottom of
# the file (under FRAGALYSIS SETTINGS).
#
# Guidance for variables: -
#
# 1.    Everything *MUST* have a default value, this file should not raise an exception
#       if a value cannot be found in the environment, that's the role of the
#       application code.
#
# 2.    The constant used to hold the environment variable *SHOULD* match the
#       environment variable's name. i.e. the "DEPLOYMENT_MODE" environment variable's
#       value *SHOULD* be found in 'settings.DEPLOYMENT_MODE' variable.
#
# 3.    In the FRAGALYSIS section, document the variable's purpose and the values
#       it can take in the comments. If there are dependencies or "gotchas"
#       (i.e. changing its value after deployment) then these should be documented.
#
# Providing run-time values for variables: -
#
# The environment variable values are set using either a 'docker-compose' file
# (when used for local development) or, more typically, via an "Ansible variable"
# provided by the "Ansible playbook" that's responsible for deploying the stack.
#
# Many (not all) of the environment variables are made available
# for deployment using an Ansible playbook variable, explained below.
#
# 1.    Ansible variables are lower-case and use "snake case".
#
# 2.    Ansible variables that map directly to environment variables in this file
#       use the same name as the environment variable and are prefixed with
#       "stack_". For example the "DEPLOYMENT_MODE" environment variable
#       can be set using the "stack_deployment_mode" variable.
#
# 3.    Variables are declared using the 'EXTRA VARIABLES' section of the corresponding
#       AWX "Job Template".
#
# IMPORTANTLY: For a description of an environment variable (setting) and its value
#              you *MUST* consult the comments in this file ("settings.py"), and *NOT*
#              the Ansible playbook. "settings.py" is the primary authority for the
#              configuration of the Fragalysis Stack.
#
# Ansible variables are declared in "roles/fragalysis-stack/defaults/main.yaml"
# or "roles/fragalysis-stack/vars/main.yaml" of the playbook repository
# https://github.com/xchem/fragalysis-stack-kubernetes
#
# For more information on "settings.py", see
# https://docs.djangoproject.com/en/3.2/topics/settings/
#
# For the full list of Django-related settings and their values, see
# https://docs.djangoproject.com/en/3.2/ref/settings/

import os
import sys
from datetime import timedelta
from typing import List, Optional

import sentry_sdk
from sentry_sdk.integrations.celery import CeleryIntegration
from sentry_sdk.integrations.django import DjangoIntegration
from sentry_sdk.integrations.excepthook import ExcepthookIntegration
from sentry_sdk.integrations.redis import RedisIntegration

# --------------------------------------------------------------------------------------
# DJANGO SETTINGS
# --------------------------------------------------------------------------------------

ALLOWED_HOSTS = ["*"]

# AnonymousUser should be the first record inserted into the auth_user table.
ANONYMOUS_USER = 1

AUTHENTICATION_BACKENDS = (
    "django.contrib.auth.backends.ModelBackend",
    "fragalysis.auth.KeycloakOIDCAuthenticationBackend",
    "guardian.backends.ObjectPermissionBackend",
)

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator"
    },
    {"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator"},
    {"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator"},
    {"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator"},
]

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# CELERY STUFF
CELERY_ACCEPT_CONTENT = ["application/json"]
CELERY_BROKER_CONNECTION_RETRY_ON_STARTUP = True
CELERY_BROKER_URL = os.environ.get("CELERY_BROKER_URL", "redis://redis:6379/")
CELERY_RESULT_BACKEND = os.environ.get("CELERY_RESULT_BACKEND", "redis://redis:6379/0")
CELERY_RESULT_BACKEND_ALWAYS_RETRY = True
CELERY_RESULT_EXPIRES = timedelta(days=15)
CELERY_TASK_ALWAYS_EAGER = os.environ.get(
    "CELERY_TASK_ALWAYS_EAGER", "False"
).lower() in ["true", "yes"]
CELERY_WORKER_HIJACK_ROOT_LOGGER = False

# SECURITY WARNING: don't run with DUBUG turned on in production!
DEBUG = os.environ.get("DEBUG_FRAGALYSIS") == "True"

DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

# Application definition
INSTALLED_APPS = [
    "django.contrib.auth",
    "django.contrib.admin",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    # 3rd
    "django_celery_beat",
    "django_celery_results",
    # My own apps
    "scoring",
    "network",
    "viewer",
    "api",
    "hypothesis",
    "hotspots",
    "media_serve",
    "service_status.apps.ServiceStatusConfig",
    # The XChem database model
    "xchem_db",
    # My utility apps
    "bootstrap3",
    "guardian",
    "graphene_django",
    "django_filters",
    "mozilla_django_oidc",  # Load after auth
    "django_extensions",
    "rest_framework",
    "rest_framework.authtoken",
    "rest_framework_swagger",
    "webpack_loader",
    "django_cleanup",
    "simple_history",
    "django_prometheus",
]

LANGUAGE_CODE = "en-us"

# Swagger logging / logout
LOGIN_URL = "/accounts/login/"
LOGOUT_URL = "/accounts/logout/"
# LOGIN_REDIRECT_URL = "<URL path to redirect to after login>"
LOGIN_REDIRECT_URL = "/viewer/react/landing"
# LOGOUT_REDIRECT_URL = "<URL path to redirect to after logout.
# Must be in keycloak call back if used>"
LOGOUT_REDIRECT_URL = "/viewer/react/landing"

MIDDLEWARE = [
    "django_prometheus.middleware.PrometheusBeforeMiddleware",
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
    "mozilla_django_oidc.middleware.SessionRefresh",
    "django_prometheus.middleware.PrometheusAfterMiddleware",
]

PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, ".."))

REST_FRAMEWORK = {
    "DEFAULT_FILTER_BACKENDS": ("django_filters.rest_framework.DjangoFilterBackend",),
    "DEFAULT_PAGINATION_CLASS": "rest_framework.pagination.LimitOffsetPagination",
    "PAGE_SIZE": 5000,
    "DEFAULT_VERSIONING_CLASS": "rest_framework.versioning.QueryParameterVersioning",
    "DEFAULT_AUTHENTICATION_CLASSES": [
        "rest_framework.authentication.SessionAuthentication",
        "mozilla_django_oidc.contrib.drf.OIDCAuthentication",
        "rest_framework.authentication.BasicAuthentication",
    ],
}

ROOT_URLCONF = "fragalysis.urls"

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = os.environ.get(
    "WEB_DJANGO_SECRET_KEY", "8flmz)c9i!o&f1-moi5-p&9ak4r9=ck$3!0y1@%34p^(6i*^_9"
)

if SENTRY_DNS := os.environ.get("FRAGALYSIS_BACKEND_SENTRY_DNS"):
    # By default only call sentry in staging/production
    sentry_sdk.init(
        dsn=SENTRY_DNS,
        integrations=[
            DjangoIntegration(),
            CeleryIntegration(),
            RedisIntegration(),
            ExcepthookIntegration(always_run=True),
        ],
        # If you wish to associate users to errors (assuming you are using
        # django.contrib.auth) you may enable sending PII data.
        send_default_pii=True,
    )

STATIC_ROOT = os.path.join(PROJECT_ROOT, "static")
STATICFILES_DIRS = [os.path.join(BASE_DIR, "fragalysis", "../viewer/static")]
STATICFILES_FINDERS = (
    "django.contrib.staticfiles.finders.FileSystemFinder",
    "django.contrib.staticfiles.finders.AppDirectoriesFinder",
)

USE_X_FORWARDED_HOST = True
SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")

# A list of identifiers of messages generated by the system check framework
# that we wish to permanently acknowledge and ignore.
# Silenced checks will not be output to the console.
#
# fields.W342   Is issued for the xchem-db package.
#               The hint is "ForeignKey(unique=True) is usually better served by a OneToOneField."
SILENCED_SYSTEM_CHECKS = [
    "fields.W342",
]

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ]
        },
    }
]

TIME_ZONE = "UTC"

# mozilla_django_oidc.
# See: https://mozilla-django-oidc.readthedocs.io/en/stable/
# Before you can configure your application, you need to set up a client with
# an OpenID Connect provider (OP). You’ll need to set up a different client for
# every environment you have for your site. For example, if your site has a -dev,
# -stage, and -prod environments, each of those has a different hostname and thus you
# need to set up a separate client for each one.
# You need to provide your OpenID Connect provider (OP) the callback url for your site.
# The URL path for the callback url is /oidc/callback/.
#
# Here are examples of callback urls:
#
#   http://127.0.0.1:8000/oidc/callback/ – for local development
#   https://myapp-dev.example.com/oidc/callback/ – -dev environment for myapp
#   https://myapp.herokuapps.com/oidc/callback/ – my app running on Heroku
#
# The OpenID Connect provider (OP) will then give you the following:
#
#   a client id (OIDC_RP_CLIENT_ID)
#   a client secret (OIDC_RP_CLIENT_SECRET)

# Keycloak mozilla_django_oidc settings (openid provider = OP).
# These should be environment variables - not checked in
OIDC_RP_CLIENT_ID = os.environ.get("OIDC_RP_CLIENT_ID", "fragalysis-local")
OIDC_RP_CLIENT_SECRET = os.environ.get("OIDC_RP_CLIENT_SECRET")
OIDC_KEYCLOAK_REALM = os.environ.get(
    "OIDC_KEYCLOAK_REALM", "https://keycloak.xchem-dev.diamond.ac.uk/auth/realms/xchem"
)

# Squonk2 Account Server and Data Manager Client IDs
OIDC_AS_CLIENT_ID: str = os.environ.get("OIDC_AS_CLIENT_ID", "")
OIDC_DM_CLIENT_ID: str = os.environ.get("OIDC_DM_CLIENT_ID", "")

# OIDC_OP_AUTHORIZATION_ENDPOINT = "<URL of the OIDC OP authorization endpoint>"
OIDC_OP_AUTHORIZATION_ENDPOINT = os.path.join(
    OIDC_KEYCLOAK_REALM, "protocol/openid-connect/auth"
)
# OIDC_OP_TOKEN_ENDPOINT = "<URL of the OIDC OP token endpoint>"
OIDC_OP_TOKEN_ENDPOINT = os.path.join(
    OIDC_KEYCLOAK_REALM, "protocol/openid-connect/token"
)
# OIDC_OP_USER_ENDPOINT = "<URL of the OIDC OP userinfo endpoint>"
OIDC_OP_USER_ENDPOINT = os.path.join(
    OIDC_KEYCLOAK_REALM, "protocol/openid-connect/userinfo"
)
# OIDC_OP_JWKS_ENDPOINT = "<URL of the OIDC OP certs endpoint>"
# This is required when using RS256.
OIDC_OP_JWKS_ENDPOINT = os.path.join(
    OIDC_KEYCLOAK_REALM, "protocol/openid-connect/certs"
)
# OIDC_OP_LOGOUT_ENDPOINT = "<URL of the OIDC OP certs endpoint>"
# This is required when using RS256.
OIDC_OP_LOGOUT_ENDPOINT = os.path.join(
    OIDC_KEYCLOAK_REALM, "protocol/openid-connect/logout"
)

# Override method to also log user out from Keycloak as well as Django.
# If desired, this should be set to "fragalysis.views.keycloak_logout"
OIDC_OP_LOGOUT_URL_METHOD = os.environ.get("OIDC_OP_LOGOUT_URL_METHOD")

# After much trial and error
# Using RS256 + JWKS Endpoint seems to work with no value for OIDC_RP_IDP_SIGN_KEY
# seems to work for authentication. Trying HS256 produces a "JWS token verification failed"
# error for some reason.
OIDC_RP_SIGN_ALGO = "RS256"
OIDC_STORE_ACCESS_TOKEN = True
OIDC_STORE_ID_TOKEN = True

# SessionRefresh configuration.
# There's only one item - the token expiry period, with a default of 15 minutes.
# The default is 15 minutes if you don't set this value.
TOKEN_EXPIRY_MINUTES = os.environ.get("OIDC_RENEW_ID_TOKEN_EXPIRY_MINUTES", "15")
OIDC_RENEW_ID_TOKEN_EXPIRY_SECONDS = int(TOKEN_EXPIRY_MINUTES) * 60

WSGI_APPLICATION = "fragalysis.wsgi.application"

DATABASE_ROUTERS = ["xchem_db.routers.AuthRouter"]

DATABASES = {
    "default": {
        "ENGINE": "django_prometheus.db.backends.postgresql",
        "NAME": os.environ.get("POSTGRESQL_DATABASE", "frag"),
        "USER": os.environ.get("POSTGRESQL_USER", "fragalysis"),
        "PASSWORD": os.environ.get("POSTGRESQL_PASSWORD", "fragalysis"),
        "HOST": os.environ.get("POSTGRESQL_HOST", "database"),
        "PORT": os.environ.get("POSTGRESQL_PORT", 5432),
    }
}

CHEMCENTRAL_DB_NAME = os.environ.get("CHEMCENT_DB_NAME", "UNKNOWN")
if CHEMCENTRAL_DB_NAME != "UNKNOWN":
    DATABASES["chemcentral"] = {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": CHEMCENTRAL_DB_NAME,
        "USER": os.environ.get("CHEMCENT_DB_USER", "postgres"),
        "PASSWORD": os.environ.get("CHEMCENT_DB_PASSWORD", "postgres"),
        "HOST": os.environ.get("CHEMCENT_DB_HOST", "postgres"),
        "PORT": 5432,
    }

USE_I18N = True
USE_L10N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
STATIC_URL = "/static/"
MEDIA_ROOT = "/code/media/"
MEDIA_URL = "/media/"

WEBPACK_LOADER = {
    "DEFAULT": {
        "BUNDLE_DIR_NAME": "bundles/",
        "STATS_FILE": os.path.join(BASE_DIR, "frontend", "webpack-stats.json"),
    }
}

GRAPHENE = {"SCHEMA": "fragalysis.schema.schema"}  # Where your Graphene schema lives

GRAPH_MODELS = {"all_applications": True, "group_models": True}

# email settings for upload key stuff
EMAIL_BACKEND = "django.core.mail.backends.smtp.EmailBackend"
if EMAIL_HOST_USER := os.environ.get("EMAIL_USER"):
    EMAIL_HOST = os.environ.get("EMAIL_HOST", "smtp.gmail.com")
    EMAIL_USE_TLS = os.environ.get("EMAIL_USE_TLS", True)
    EMAIL_PORT = os.environ.get("EMAIL_PORT", 587)
    EMAIL_HOST_PASSWORD = os.environ.get("EMAIL_PASSWORD")

# Configure django logging.
# We provide a standard formatter that emits a timestamp, the module issuing the log
# and the level name, a little like this...
#
#   2022-05-16T09:04:29 django.request ERROR # Internal Server Error: /viewer/react/landing
#
# We provide a console and rotating file handler
# (50Mi of logging in 10 files of 5M each),
# with the rotating file handler typically used for everything.
DISABLE_LOGGING_FRAMEWORK = os.environ.get(
    "DISABLE_LOGGING_FRAMEWORK", "no"
).lower() in ["yes"]
LOGGING_FRAMEWORK_ROOT_LEVEL = os.environ.get("LOGGING_FRAMEWORK_ROOT_LEVEL", "DEBUG")
if not DISABLE_LOGGING_FRAMEWORK:
    LOGGING = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "simple": {
                "format": "%(asctime)s %(name)s.%(funcName)s():%(lineno)s %(levelname)s # %(message)s",
                "datefmt": "%Y-%m-%dT%H:%M:%S%z",
            }
        },
        "handlers": {
            "console": {
                "level": "DEBUG",
                "class": "logging.StreamHandler",
                "stream": sys.stdout,
                "formatter": "simple",
            },
            "rotating": {
                "level": "DEBUG",
                "class": "logging.handlers.RotatingFileHandler",
                "maxBytes": 5_000_000,
                "backupCount": 10,
                "filename": os.path.join(BASE_DIR, "logs/backend.log"),
                "formatter": "simple",
            },
            'service_status': {
                'level': 'DEBUG',
                'class': 'logging.handlers.RotatingFileHandler',
                'filename': os.path.join(BASE_DIR, "logs/service_status.log"),
                'formatter': 'simple',
                "maxBytes": 5_000_000,
                "backupCount": 10,
            },
        },
        "root": {
            "level": LOGGING_FRAMEWORK_ROOT_LEVEL,
            "handlers": ["console", "rotating"],
        },
        'loggers': {
            'api.security': {'level': 'INFO'},
            'asyncio': {'level': 'WARNING'},
            'celery': {'level': 'INFO'},
            'django': {'level': 'ERROR'},
            'mozilla_django_oidc': {'level': 'WARNING'},
            'urllib3': {'level': 'WARNING'},
            'paramiko': {'level': 'WARNING'},
            'service_status': {
                'handlers': ['service_status', 'console'],
                'level': 'DEBUG',
                'propagate': False,
            },
        },
    }

# --------------------------------------------------------------------------------------
# FRAGALYSIS SETTINGS
# --------------------------------------------------------------------------------------
# With comprehensive comments where necessary to explain the setting's values.

# The deployment mode.
# Controls the behaviour of the application (it's strictness to errors etc).
# Typically one of "DEVELOPMENT" or "PRODUCTION".
# see api.utils for the 'deployment_mode_is_production()' function.
DEPLOYMENT_MODE: str = os.environ.get("DEPLOYMENT_MODE", "production").upper()

# Authentication check when uploading files.
# This can be switched off to simplify development testing if required.
# It's asserted as True for 'production' mode.
AUTHENTICATE_UPLOAD: bool = True
if os.environ.get("AUTHENTICATE_UPLOAD") == "False":
    assert DEPLOYMENT_MODE != "PRODUCTION"
    AUTHENTICATE_UPLOAD = False

COMPUTED_SET_MEDIA_DIRECTORY: str = "computed_set_data"

# The following (part of m2ms-1385) is used to prevent the
# 'restrict-to-membership' check in security.py - something that is designed to prevent
# uploading to public proposals unless the user is explicitly part of the proposal
# (according to ISPyB). This variable is used to defeat this test for situations
# when ISPyB is unavailable. It is not permitted when the DEPLOYMENT_MODE
# is 'PRODUCTION
DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP: bool = False
if os.environ.get("DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP") == "True":
    assert DEPLOYMENT_MODE != "PRODUCTION"
    DISABLE_RESTRICT_PROPOSALS_TO_MEMBERSHIP = True

# Discourse settings for API calls to Discourse Platform.
DISCOURSE_PARENT_CATEGORY: str = "Fragalysis targets"
DISCOURSE_USER: str = "fragalysis"
DISCOURSE_HOST: str = os.environ.get("DISCOURSE_HOST", "")
# If a DISCOURSE_API_KEY is not set the backend will assume
# that a Discourse server is not available.
# Note that this can be obtained from discourse for the desired environment.
DISCOURSE_API_KEY: str = os.environ.get("DISCOURSE_API_KEY", "")
# This suffix can be set to that the different development environments posting
# to the same Discourse server can "automatically" generate different category/post
# titles - hopefully reducing confusion. It will be appended at category or post-title,
# e.g. "Mpro-duncan", "Mpro-staging" etc. Note that it is for dev systems.
# It is not required on production because production will have a
# dedicated Discourse server.
DISCOURSE_DEV_POST_SUFFIX: str = os.environ.get("DISCOURSE_DEV_POST_SUFFIX", "")

# Some Squonk2 developer/debug variables.
# Unused in production.
DUMMY_TARGET_TITLE: str = os.environ.get("DUMMY_TARGET_TITLE", "")
DUMMY_USER: str = os.environ.get("DUMMY_USER", "")
DUMMY_TAS: str = os.environ.get("DUMMY_TAS", "")

# Do we enable the collection and presentation
# of the availability of underlying services?
# A colon (:) separated list of services to enable.
# See "viewer/services.py" for the full list of supported services.
ENABLE_SERVICE_STATUS: str = os.environ.get("ENABLE_SERVICE_STATUS", "")

# What infection have been set?
# "Infections" are  built-in faults that can be induced by providing their names.
# Typically these are "hard to reproduce" errors that are useful for testing.
# The names are provided in a comma-separated list in this variable.
# The full set of supported names can be used can be found in "api/infections.py"
INFECTIONS: str = os.environ.get("INFECTIONS", "").lower()

# The ISpyB database settings.
# Can be used in conjunction with SSH settings (later in this file)
ISPYB_USER: str = os.environ.get("ISPYB_USER", "")
ISPYB_PASSWORD: str = os.environ.get("ISPYB_PASSWORD", "")
ISPYB_HOST: str = os.environ.get("ISPYB_HOST", "")
ISPYB_PORT: str = os.environ.get("ISPYB_PORT", "")

# An optional URL that identifies the URL to a prior stack.
# If set, it's typically something like "https://fragalysis.diamond.ac.uk".
# It can be blank, indicating there is no legacy service.
LEGACY_URL: str = os.environ.get("LEGACY_URL", "")

NEOMODEL_NEO4J_BOLT_URL: str = os.environ.get(
    "NEO4J_BOLT_URL", "bolt://neo4j:test@neo4j:7687"
)

# The graph (neo4j) database settings.
# The query provides the graph endpoint, typically a service in a kubernetes namespace
# like 'graph.graph-a.svc' and the 'auth' provides the graph username and password.
NEO4J_QUERY: str = os.environ.get("NEO4J_QUERY", "neo4j")
NEO4J_AUTH: str = os.environ.get("NEO4J_AUTH", "neo4j/neo4j")

# Does it look like we're running in Kubernetes?
# If so, let's get the namespace we're in - it will provide
# useful discrimination material in log/metrics messages.
# If there is no apparent namespace the variable will be 'None'.
OUR_KUBERNETES_NAMESPACE: Optional[str] = None
_NS_FILENAME: str = '/var/run/secrets/kubernetes.io/serviceaccount/namespace'
if os.path.isfile(_NS_FILENAME):
    with open(_NS_FILENAME, 'rt', encoding='utf8') as ns_file:
        OUR_KUBERNETES_NAMESPACE = ns_file.read().strip()

# These flags are used in the upload_tset form as follows.
# Proposal Supported | Proposal Required | Proposal / View fields
# Y                  | Y                 | Shown / Required
# Y                  | N                 | Shown / Optional
# N                  | N                 | Not Shown
PROPOSAL_SUPPORTED: bool = True
PROPOSAL_REQUIRED: bool = True

# Are any public target access strings defined?
# If so they'll be in the PUBLIC_TAS variable as a comma separated list.
PUBLIC_TAS: str = os.environ.get("PUBLIC_TAS", "")
PUBLIC_TAS_LIST: List[str] = PUBLIC_TAS.split(",") if PUBLIC_TAS else []

# Security/access control connector.
# Currently one of 'ispyb' or 'ssh_ispyb'.
SECURITY_CONNECTOR: str = os.environ.get("SECURITY_CONNECTOR", "ispyb").lower()
# Number of minutes to cache security information for a user.
# Set to '0' to disable caching.
SECURITY_CONNECTOR_CACHE_MINUTES: int = int(
    os.environ.get("SECURITY_CONNECTOR_CACHE_MINUTES", "2")
)

# An SSH host.
# Used in the security module in conjunction with ISPyB settings.
# The SSH_PRIVATE_KEY_FILENAME value will be used if there is no SSH_PASSWORD.
SSH_HOST: str = os.environ.get("SSH_HOST", "")
SSH_USER: str = os.environ.get("SSH_USER", "")
SSH_PASSWORD: str = os.environ.get("SSH_PASSWORD", "")
SSH_PRIVATE_KEY_FILENAME: str = os.environ.get("SSH_PRIVATE_KEY_FILENAME", "")

# The maximum length of the 'slug' used for names this Fragalysis will create.
#
# Squonk2 variables are generally used by the 'squonk2_agent.py' module
# in the 'viewer' package.
SQUONK2_MAX_SLUG_LENGTH: int = 10

# Where the Squonk2 logic places its files in Job containers.
SQUONK2_MEDIA_DIRECTORY: str = "fragalysis-files"
# The Squonk2 DataManger UI endpoint to obtain Job Instance information.
SQUONK2_INSTANCE_API: str = "data-manager-ui/results/instance/"

# The URL for the Squonk2 Account Server API.
SQUONK2_ASAPI_URL: str = os.environ.get("SQUONK2_ASAPI_URL", "")
# The URL for the Squonk2 Data Manaqer API.
SQUONK2_DMAPI_URL: str = os.environ.get("SQUONK2_DMAPI_URL", "")
# The URL for the Squonk2 User Interface.
SQUONK2_UI_URL: str = os.environ.get("SQUONK2_UI_URL", "")
# The pre-assigned Squonk2 Account Server Organisation for the stack.
# This is created by an administrator of the Squonk2 service.
SQUONK2_ORG_UUID: str = os.environ.get("SQUONK2_ORG_UUID", "")
# The Account Server Unit billing day 9for all products (projects) that are created.
# It's a day of the month (1..27).
SQUONK2_UNIT_BILLING_DAY: str = os.environ.get("SQUONK2_UNIT_BILLING_DAY", "")
# The Squonk2 Account Server product "flavour" created for Jobs (products/projects).
# It's usually one of "GOLD", "SILVER" or "BRONZE".
SQUONK2_PRODUCT_FLAVOUR: str = os.environ.get("SQUONK2_PRODUCT_FLAVOUR", "")
# A short slug used when creating Squonk2 objects for this stack.
# This must be unique across all stacks that share the same Squonk2 service.
SQUONK2_SLUG: str = os.environ.get("SQUONK2_SLUG", "")[:SQUONK2_MAX_SLUG_LENGTH]
# The pre-assigned Squonk2 Account Server Organisation owner and password.
# This account is used to create Squonk2 objects for the stack.
SQUONK2_ORG_OWNER: str = os.environ.get("SQUONK2_ORG_OWNER", "")
SQUONK2_ORG_OWNER_PASSWORD: str = os.environ.get("SQUONK2_ORG_OWNER_PASSWORD", "")
# Do we verify Squonk2 SSL certificates ("yes" or "no").
SQUONK2_VERIFY_CERTIFICATES: str = os.environ.get("SQUONK2_VERIFY_CERTIFICATES", "")

TARGET_LOADER_MEDIA_DIRECTORY: str = "target_loader_data"

# A warning messages issued by the f/e.
# Used, if set, to populate the 'target_warning_message' context variable
TARGET_WARNING_MESSAGE: str = os.environ.get("TARGET_WARNING_MESSAGE", "")

# The Target Access String (TAS) Python regular expression.
# The Project title (the TAS) must match this expression to be valid.
# See api/utils.py validate_tas() for the current implementation.
# To simplify error messages when the match fails you can also
# add an error message.
TAS_REGEX: str = os.environ.get("TAS_REGEX", r"^(lb\d{5})(-(\d+)){0,1}$")
TAS_REGEX_ERROR_MSG: str = os.environ.get(
    "TAS_REGEX_ERROR_MSG",
    "Must begin 'lb' followed by 5 digits, optionally followed by a hyphen and a number.",
)

# Version variables.
# These are set by the Dockerfile in the fragalysis-stack repository
# and controlled by the CI process, i.e. they're not normally set by a a user.
BE_NAMESPACE: str = os.environ.get("BE_NAMESPACE", "undefined")
BE_IMAGE_TAG: str = os.environ.get("BE_IMAGE_TAG", "undefined")
FE_NAMESPACE: str = os.environ.get("FE_NAMESPACE", "undefined")
FE_IMAGE_TAG: str = os.environ.get("FE_IMAGE_TAG", "undefined")
STACK_NAMESPACE: str = os.environ.get("STACK_NAMESPACE", "undefined")
STACK_VERSION: str = os.environ.get("STACK_VERSION", "undefined")
