# Integration tests — developer guide

This guide covers the **API integration tests**: tests that drive the real DRF
API end-to-end against a *running* stack (a real database, a redis broker and a
celery worker), using large test archives pulled from a public-read S3 bucket.

It complements the fast **unit tests** (`viewer/tests/`, run with a bare
`pytest` against a single Postgres container — see the project `README` /
`CLAUDE.md`). The two layers are deliberately separate:

| | Unit tests | Integration tests |
|---|---|---|
| Runner | host `pytest` (`run-unit-tests.sh`) | `pytest` inside the backend image |
| Needs | one Postgres container | full stack: db + redis + **worker** + backend |
| Celery | eager (in-process) | **real** broker + worker (async) |
| Data | hand-built fixtures | large archive from S3, loaded via the API |
| Speed | seconds | minutes (download + load) |
| Default `pytest`? | yes | **no** — deselected and env-gated |

The async upload flow hands work to a Celery worker over the broker, so it is
only exercised faithfully against a real worker — *not* under eager mode. That
is why there is no in-process counterpart to the integration test.

---

## The external-data mechanism

Realistic target archives are far too large to commit, so they live in a
**public-read** S3 bucket (anonymous `GET` only; no `ListBucket`) and are
downloaded over plain HTTPS at test time.

Two environment variables gate and select the data:

| Variable | Meaning | Example |
|---|---|---|
| `UNIT_TEST_BUCKET_AND_PATH` | `s3://bucket/prefix` root of the test data | `s3://im-fragalysis-backend/unit-test` |
| `UNIT_TEST_DATA_IDENTIFIER` | which manifest entry to use | `ALPHA` |

When either is unset, the tests **skip** (so the default `poetry run pytest` is
unaffected). The download is anonymous, so **no AWS credentials are needed.**

### Files

| Path | Role |
|---|---|
| `viewer/tests/external_data.py` | Shared helpers: the `requires_external_data` skip marker, `s3://`→`https` conversion, a streamed anonymous `urllib3` download (raises on non-200), and the manifest loader. |
| `viewer/tests/test_data/api/<endpoint>/manifest.yaml` | Per-endpoint manifest: for each identifier, the object(s) to upload and the results to expect. |
| `viewer/tests/test_external_data.py` | Unconditional unit tests for the pure helper logic (URL parsing, env gating). These run in the normal suite. |
| `viewer/tests/test_api_upload_target_experiments_async.py` | The integration test itself (marked `integration`). |
| `docker-compose.integration.yml` | The self-contained stack the test runs against. |
| `.github/workflows/integration-tests.yaml` | Reusable CI workflow that stands up the stack and runs the test. |

### S3 layout

The bucket uses a **flat** layout — the archive sits directly under the
identifier, with no per-TAS sub-directory:

```
<UNIT_TEST_BUCKET_AND_PATH>/api/<endpoint>/<IDENTIFIER>/<file>
e.g.  s3://im-fragalysis-backend/unit-test/api/upload_target_experiments/ALPHA/A71EV2A-20260603.tgz
```

### The manifest

`viewer/tests/test_data/api/upload_target_experiments/manifest.yaml`:

```yaml
ALPHA:
  uploads:
    - tas: lb32627-66            # target_access_string, passed verbatim to the API
      file: A71EV2A-20260603.tgz # object under the identifier directory
  expect:
    target_title: A71EV2A
    targets: 1
    experiments: 7642
    site_observations: 67
    poses: 62
    target_experiment_uploads: 1
```

`target_title` and `targets` are always asserted. Each remaining count is only
asserted when **non-null**, so a new identifier can start with `null`
placeholders and gain assertions as the real numbers are filled in.

---

## How the async test works

`test_api_upload_target_experiments_async.py` is marked `@pytest.mark.integration`
**and** gated on `INTEGRATION_BASE_URL` (the base URL of the running stack). It
is driven purely over HTTP with `urllib3` — it owns no Django test DB; the
running stack owns the one real database. For each manifest upload it:

1. **Downloads** the archive from S3 (anonymous).
2. **POSTs** it to `/api/upload_target_experiments/` (multipart). The stack runs
   with `AUTHENTICATE_UPLOAD=False`, so the upload is accepted without a
   pre-existing proposal and returns **202** with a `task_status_url`.
3. **Polls** `task_status_url` until the task reports `finished`, then asserts
   `status == "SUCCESS"`. The poll tolerates transient non-200s: while the load
   runs, `task_status` briefly returns **404** ("Proposal not found") because
   the proposal's `Project` is not committed until partway through the load.
4. **Asserts** the GET endpoints against the manifest `expect`:
   `/api/targets/`, `/api/experiments/`, `/api/site_observations/`,
   `/api/poses/`, `/api/target_experiment_uploads/`.

Visibility for the anonymous GET/poll comes from the stack's `PUBLIC_TAS`, which
publishes the loaded proposal (the loaded `Project.title` equals the TAS) so the
data is visible without membership.

### pytest wiring

`pyproject.toml` (`[tool.pytest.ini_options]`) registers the marker and
deselects it by default:

```toml
addopts = "--reuse-db --strict-markers -m \"not integration\" --cov=viewer ..."
markers = ["integration: needs the live container stack (redis + celery worker)"]
```

So neither the host run nor the standard CI `pytest` job ever collects the
integration test — it runs **only** inside the integration stack, which passes
`-m integration`.

---

## The stack (`docker-compose.integration.yml`)

| Service | What it is |
|---|---|
| `database` | Postgres (pgvector). Uses an **ephemeral named volume**, so every run starts clean. |
| `redis` | The Celery broker + result backend (ephemeral volume). |
| `backend` | gunicorn + nginx on `:80`. `CELERY_TASK_ALWAYS_EAGER=False`, `AUTHENTICATE_UPLOAD=False`, `PUBLIC_TAS=lb32627-66`, broker/result backend pointed at `redis`. |
| `celery_worker` | A **real** worker (`celery -A fragalysis worker`) that processes the upload asynchronously. |
| `test` | The test runner: runs `pytest` against `INTEGRATION_BASE_URL=http://backend:80`. |

A few non-obvious points baked into the compose file:

- **Shared media volume.** The upload view writes the bundle to `MEDIA_ROOT`
  (`/code/media/tmp/...`) and the *worker* (a separate container) reads it back,
  so `backend` and `celery_worker` share a `media` volume.
- **Anonymous committer.** With no authenticated user, the loader falls back to
  `settings.ANONYMOUS_USER` (pk 1), which is the `admin` superuser that
  `launch-stack.sh` creates on first boot.
- **`DJANGO_SETTINGS_MODULE=fragalysis.settings` on the `test` service.**
  `pyproject` points pytest at `tests.test_settings`, but the top-level `tests/`
  package is `.dockerignore`'d out of the image. The integration test is pure
  HTTP, but pytest still imports `viewer/tests/conftest.py` (which imports Django
  models), so Django must be configured — using the app's own settings module,
  which *is* present in the image.
- **Test toolchain installed at runtime.** The production image carries only the
  `main` dependencies, so the `test` service `pip install`s `pytest`,
  `pytest-django` and `pytest-cov` before running.

---

## Running locally

You need a backend **image built from your checkout** (the `test` service runs
pytest from inside the image, so it must contain your test files). Build it with
the normal build-stack flow, or directly:

```bash
docker build -t xchem/fragalysis-backend:dev .
```

Then bring the stack up and let the `test` service drive it:

```bash
BE_NAMESPACE=xchem BE_IMAGE_TAG=dev \
UNIT_TEST_BUCKET_AND_PATH=s3://im-fragalysis-backend/unit-test \
UNIT_TEST_DATA_IDENTIFIER=ALPHA \
docker compose -f docker-compose.integration.yml up \
  --abort-on-container-exit --exit-code-from test
```

`--exit-code-from test` makes the whole `up` exit with the test's status. Tear
down (and wipe the ephemeral volumes) with:

```bash
docker compose -f docker-compose.integration.yml down -v
```

### Apple-Silicon gotcha (pgvector SIGILL)

On arm64 (Apple Silicon) the database image's HNSW index on
`viewer_atomcoordinates.coords` (`pgvector_coord_index`) crashes Postgres with
*"Illegal instruction"* during per-row index maintenance — the load aborts and
the DB enters recovery. This is a CPU/SIMD mismatch in that image's pgvector
build; **CI runs on linux/amd64 and is unaffected.** The normal unit suite never
hits it because it doesn't bulk-insert coordinate vectors.

To load locally on Apple Silicon, drop that index after the backend has migrated
but before the load (row counts don't depend on it):

```bash
docker compose -f docker-compose.integration.yml up -d --wait database redis backend celery_worker
docker exec database psql -U postgres -d frag -c "DROP INDEX IF EXISTS pgvector_coord_index;"
docker compose -f docker-compose.integration.yml run --rm --no-deps test
```

---

## Running on CI

The stack is wrapped in a **reusable** workflow,
`.github/workflows/integration-tests.yaml` (`on: workflow_call`):

- **Inputs:** `be_namespace`, `be_image_tag`, `unit_test_bucket_and_path`,
  `unit_test_data_identifier`. **No secrets** (the download is anonymous).
- **Steps:** check out → `docker compose ... up --abort-on-container-exit
  --exit-code-from test` → always dump logs → `down -v`.

It is wired into all three build workflows (`build-dev`, `build-staging`,
`build-production`) as a downstream job:

```yaml
integration-tests:
  needs: build
  if: needs.build.outputs.push == 'true'   # the test stack pulls the image by tag
  uses: ./.github/workflows/integration-tests.yaml
  with:
    be_namespace: ${{ needs.build.outputs.be_namespace }}
    be_image_tag: ${{ needs.build.outputs.tag }}
    unit_test_bucket_and_path: ${{ needs.build.outputs.unit_test_bucket_and_path }}
    unit_test_data_identifier: ${{ needs.build.outputs.unit_test_data_identifier }}
```

The job runs against the **freshly built image** and only when that image was
pushed (otherwise the tag the stack pulls would not exist). The bucket and
identifier are defined once as workflow-level `env` (`UNIT_TEST_BUCKET_AND_PATH`,
`UNIT_TEST_DATA_IDENTIFIER`) and surfaced as build-job outputs.

---

## Adding a new dataset

1. Upload the archive to the bucket at
   `api/<endpoint>/<IDENTIFIER>/<file>` (flat — no per-TAS directory).
2. Add an entry to the endpoint's `manifest.yaml` with the `tas` and `file`, and
   the `expect` results. Start the counts at `null` if you don't know them yet.
3. Capture the real `expect` counts by running the stack once and reading the
   endpoint `count`s back (this is how the `ALPHA` numbers were obtained) — then
   replace the `null`s. A non-null count turns on that endpoint's assertion.
4. Point CI/local runs at it with `UNIT_TEST_DATA_IDENTIFIER=<IDENTIFIER>`.

---

## Troubleshooting

- **`task_status` 404 during the load** — expected and tolerated; the proposal's
  `Project` is not committed until partway through. The poll keeps trying until
  the deadline.
- **A re-run is rejected as a duplicate upload** ("next version should be 2") —
  the database was not reset between runs. The stack uses ephemeral named
  volumes; make sure you `down -v` (the CI workflow always does).
- **`No module named 'tests'`** in the `test` service — `DJANGO_SETTINGS_MODULE`
  must be `fragalysis.settings` (the top-level `tests/` package is not in the
  image); the compose file sets this.
- **Postgres "Illegal instruction" / recovery mode on arm64** — the pgvector
  HNSW index; see the Apple-Silicon gotcha above.
