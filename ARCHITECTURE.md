# Architecture

This document describes the *big-picture* design of the fragalysis-backend — the
parts that span multiple files and are not obvious from any single module. For
setup/commands see `README.md`; for working conventions see `CLAUDE.md`. Historical
design rationale for individual features lives as PDFs/PlantUML in `design_docs/`.

## The stack and where this fits

The backend is one component of the deployed **Fragalysis stack**. At runtime it
talks to four external services:

- **PostgreSQL** — the primary relational store (with `pgvector` for embeddings).
- **Neo4j** — a graph database, used for fragment-network queries (`network/`).
- **Redis** — Celery broker/result backend and Django cache.
- **Squonk2** (Data Manager + Account Server) — an external compute platform that
  runs computational jobs; reached over HTTP via the `squonk2` Python client.

The process model inside the container (`launch-stack.sh`): migrations →
load fixtures → `collectstatic` → create superuser → **gunicorn** (the Django app)
behind **nginx**. On the lead pod (`stack-0`) it also starts the service health-check
scheduler and the download-cleanup scheduler. Celery **workers** and **beat** run as
separate processes (`launch-worker.sh`, `launch-beat.sh`).

## Django project layout

`fragalysis/` is the Django *project* (settings, root URL conf, celery app, OIDC
auth). The functionality is split across several *apps* (`INSTALLED_APPS` in
`fragalysis/settings.py`), but the centre of gravity is overwhelmingly the **`viewer`**
app — `viewer/views.py`, `viewer/models.py`, `viewer/serializers.py` and
`viewer/target_loader.py` are each thousands of lines and contain most of the domain
logic. The other apps (`scoring`, `hotspots`, `hypothesis`, `network`, `media_serve`,
`service_status`) are comparatively small and focused.

URL routing starts in `fragalysis/urls.py`, which includes each app's `urls.py`. Most
data endpoints are DRF `ViewSet`s registered under `/api/`.

## Domain model

The data model (`viewer/models.py`) represents crystallographic fragment-screening
data. The core hierarchy, roughly top-down:

```
Project (a.k.a. proposal/visit — the unit of ACCESS CONTROL)
└── Target (a protein target)
    ├── ExperimentUpload (one ingestion of data; created by the target loader)
    │   └── Experiment (a crystal dataset; links Compounds + an Xtalform)
    │       └── SiteObservation  ← the central "thing a user looks at":
    │              a single observed ligand binding event
    ├── Xtalform / XtalformSite / QuatAssembly  (crystal-form description)
    ├── CanonSite / CanonSiteConf  (canonical binding sites + conformations)
    └── Pose  (groups SiteObservations of one Compound at one CanonSite)

Compound  — a unique 2D molecule (shared across targets)
```

Key relationships to know before touching models:

- **`SiteObservation`** is the hub. It points to its `Experiment`, `Compound`
  (`cmpd`), `XtalformSite`, `CanonSiteConf`, and `Pose`. Most user-facing data and
  files hang off it.
- **`CanonSite` ↔ `CanonSiteConf` ↔ `SiteObservation`** form a mutually-referential
  cluster (a canon site has a reference conf; a conf has a reference observation).
  This is why loading happens in carefully ordered passes (see below).
- **`Versionable`** is an abstract base (`CanonSite`, `CanonSiteConf`,
  `XtalformSite`, `SiteObservation` inherit it) supporting versioned/upgraded records.
- Custom **managers** live in `viewer/managers.py` (one per major model) and encode
  query logic — prefer extending these over ad-hoc querysets.
- `django-simple-history` (`HistoricalRecords`) tracks changes on several models.

**Computed sets** are a parallel, user-contributed branch: a `ComputedSet` of
`ComputedMolecule`s (3D poses produced computationally) attached to a `Target`, each
referencing a `Compound` and its inspiration `SiteObservation`s, with scores in
`ScoreDescription`/`NumericalScoreValues`/`TextScoreValues`.

**Sessions/snapshots** (`SessionProject`, `Snapshot`, `SessionActions`) persist a
user's view state of a target page. **Tags** (`Tag` subclasses in `viewer/tags.py`)
group observations and sessions.

## Access control (read this before adding any data-exposing view)

Authorisation is **proposal-based**, not row-level Django permissions. Every piece of
target data ultimately belongs to a `Project` (a "proposal"/"visit"), and a user may
only see data for proposals they are a member of, plus any proposal marked
`open_to_public` or listed in `PUBLIC_TAS`.

The mechanism is `api.security.ISPyBSafeQuerySet` (a DRF `ReadOnlyModelViewSet`
subclass):

- A view sets `filter_permissions` to the ORM path from its model to the owning
  `Project` (e.g. `"target__project"`). `get_queryset()` builds a `Q` filter limiting
  results to the user's proposals OR public proposals.
- Membership comes from one of two sources, chosen by `settings.TA_AUTH_SERVICE`:
  the external **Target-Access (TA) auth connector** (`api/ta_auth_connector.py`) when
  configured, otherwise the local Django `Project.user` relation.
- `restrict_public_to_membership` distinguishes *viewing* public data (everyone) from
  *modifying/uploading* to a public proposal (explicit members only). Use the
  `user_is_member_of_*` helpers for write paths.
- Because `ISPyBSafeQuerySet` is read-only, views needing write methods add DRF
  mixins explicitly.

When adding an endpoint that returns target-scoped data, subclass `ISPyBSafeQuerySet`
and set `filter_permissions` — do not write a plain `ModelViewSet`, or you will leak
data across proposals.

Authentication itself is **OIDC via Keycloak** (`mozilla_django_oidc`), configured in
`settings.py` and `fragalysis/auth.py`.

## Asynchronous work (Celery)

Long-running operations are Celery tasks (`viewer/tasks.py`, broker = Redis). The
major ones:

- **`load_target`** — ingest a new target data package (delegates to
  `viewer/target_loader.py`).
- **`create_download`** — build a downloadable structures archive
  (`viewer/download_structures.py`); cleaned up later by a scheduler.
- **Computed-set upload** validation/processing (`viewer/cset_upload.py`).
- **Squonk2 job** file transfer / upload / request handling
  (`viewer/squonk_job_*.py`).

Locally, tasks run **eagerly** (synchronously) unless the celery compose file is
used — see `CLAUDE.md`. The pattern in views is: kick off a task, return a task id,
and let the client poll a status endpoint.

## The target-loading pipeline

`viewer/target_loader.py` (~4k lines) is the most complex subsystem. It takes an
uploaded archive describing a target (YAML metadata + structure files + an SDF/SQLite
of observations), and populates the model hierarchy above inside a transaction. Two
things to understand:

1. **Ordered, multi-pass creation.** Because of the mutually-referential
   CanonSite/Conf/SiteObservation relationships, records are created in dependency
   order and back-references are filled in later passes.
2. **Versioned re-uploads.** Re-uploading a target produces a new `ExperimentUpload`
   and upgrades `Versionable` records rather than blindly duplicating.

After a successful load it clears the view cache (`viewer/cache.py`) and derives tags.

## Squonk2 integration

External computational jobs run on the Squonk2 platform. `viewer/squonk2_agent.py`
(`Squonk2Agent`, a process-wide singleton) wraps the Squonk2 client and maps
Fragalysis `Project`s to Squonk2 `Org`/`Unit`/`Product`/`Project` records
(`Squonk2Org`/`Squonk2Unit`/`Squonk2Project` models). The `squonk_job_*.py` modules
orchestrate transferring input files to Squonk2, requesting a job (`JobRequest`),
and handling success callbacks that ingest results back as a `ComputedSet`.

## Caching, logging, configuration

- **Cache:** per-view caching keyed appropriately; toggled by `CACHE_ENABLED` and
  cleared via `viewer/cache.py` / the `clear_cache` management command.
- **Logging:** configured in `settings.py` (console + rotating file at
  `/code/logs/backend.log`); custom adapters in `viewer/logger_adapters.py`.
- **Configuration** is environment-variable driven through `settings.py`; the file's
  header comment is the authoritative style guide. `DEPLOYMENT_MODE` and `INFECTIONS`
  (`api/infections.py`) alter behaviour for production-hardening and error-path
  testing respectively.

## Management commands

Operational/maintenance tasks live in `viewer/management/commands/` (documented in
that directory's `README.md`) and `service_status/` — e.g. `load_target`,
`clear_cache`, `start_download_cleanup`, `start_service_queries`, and various
one-off data-correction commands.
