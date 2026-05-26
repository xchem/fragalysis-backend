## The project

The Django + Django REST Framework backend for the **Fragalysis** stack — a web
application for visualising and analysing X-ray crystallographic fragment-screening
data. This repo is one of several that make up the deployed "stack" (see README's
*Background* section); it owns the database models, the REST API, data loaders, and
the Celery tasks. It is deployed as a container image to Kubernetes.

For the big-picture design (domain model, access control, async tasks, the
target-loading pipeline, Squonk2 integration) read **`ARCHITECTURE.md`** — that is
the document to consult before changing models, security, or the loader.

## Commands

Dependencies are managed with **Poetry** (`pyproject.toml` / `poetry.lock`). The
application is intended to run inside its Docker container, not on the host.

### Run locally
```bash
docker-compose up -d          # builds/starts backend + postgres + neo4j + redis
docker-compose exec backend bash   # shell into the running backend
docker-compose down
```
The API is then at `http://localhost:8080/api/`. The repo is bind-mounted into the
container at `/code`, so host edits are live (no rebuild needed). Logs:
`./data/logs/backend.log`. Locally, Celery runs **synchronously**
(`CELERY_TASK_ALWAYS_EAGER=True`); to exercise real async tasks add the celery
compose file: `docker compose -f docker-compose.yml -f docker-compose.celery.yml up`.

### Migrations
Use the dedicated migration compose file and run `makemigrations` *inside* the
container (migrations are written back to the host via the bind mount, then commit
them):
```bash
docker compose -f docker-compose-migrate.yml up -d
docker compose -f docker-compose-migrate.yml exec backend bash
python manage.py makemigrations viewer --name "descriptive_name"
```

### Lint / format
Enforced via **pre-commit**:
```bash
pre-commit run --all-files
```
Important: linting is intentionally **limited to the `viewer` app** (see
`.pre-commit-config.yaml`). Migrations are excluded from all hooks.

## Conventions and gotchas

- **Settings are environment-driven.** Almost all runtime configuration lives in
  `fragalysis/settings.py`, read from environment variables. Read the comment block
  at the top of that file for the style guide before adding a new variable. Not all
  settings are dynamic (e.g. `ALLOWED_HOSTS` is static).
- **Deployment mode** (`DEPLOYMENT_MODE` = `DEVELOPMENT` | `PRODUCTION`) changes
  behaviour — production is stricter. Use `api.utils.deployment_mode_is_production()`
  rather than reading the env var directly.
- **Access control is not optional.** API views that expose project data subclass
  `api.security.ISPyBSafeQuerySet` and declare `filter_permissions` so results are
  filtered to the user's accessible proposals. Don't bypass this when adding views —
  see `ARCHITECTURE.md` § Access control.
- **"Infections"** (`api/infections.py`) let you inject specific error paths for
  testing via the `INFECTIONS` env var; they are ignored in production mode.
- **Conventional Commits** are used for commit messages.

## Key directories

- `fragalysis/` — the Django *project* (settings, root URLs, celery, wsgi, auth).
- `viewer/` — the core app: models, DRF views/serializers, the target loader,
  computed-set upload, Squonk2 job handling, tasks, tags. The bulk of the logic.
- `api/` — cross-cutting concerns: security/access-control, infections, the target-
  access (TA) auth connector.
- `scoring/`, `hotspots/`, `hypothesis/`, `network/`, `media_serve/` — smaller apps.
- `service_status/` — scheduled external-service health checks.
- `design_docs/` — historical design documents (PDF + PlantUML) explaining *why*
  major features were built; useful background, not kept in sync with code.
- `docs/` — Sphinx source for the ReadTheDocs site.
