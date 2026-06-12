---
name: deploy-stack
description: >-
  Build a developer Fragalysis Stack (backend + stack container images)
  and then deploy it to Kubernetes by triggering the developer's AWX Job Template
---

This skill takes a developer's Fragalysis Stack from source to a running deployment
in two phases: first we **build** the backend and stack container images, then we
**deploy** the stack to Kubernetes by launching the developer's AWX Job Template.

# Gather what we need from the developer

Collect every value we need from the developer **up front**, before building
anything, so the deploy at the end never stalls waiting for credentials. There are
three things to ask for; the rest are derived.

- `STACK_NAMESPACE` — the developer's DockerHub/GitHub username, e.g.
  `alanbchristie` (the namespace the stack image is pushed to)
- `CONTROLLER_USERNAME` — the developer's AWX username, e.g. `alan`
- `CONTROLLER_PASSWORD` — the developer's AWX password

The AWX base URL is always `https://awx.xchem-dev.diamond.ac.uk`, so we don't ask
for it. Export everything once — the `docker` and `awx` CLIs read these variables
automatically: -

```
export STACK_NAMESPACE=<the developer's DockerHub/GitHub username>
export CONTROLLER_HOST=https://awx.xchem-dev.diamond.ac.uk
export CONTROLLER_USERNAME=<the developer's AWX username>
export CONTROLLER_PASSWORD=<the developer's AWX password>
```

>   Do not print `CONTROLLER_PASSWORD` back to the developer, write it to a file,
    or include it in any command output. Prefer reading it interactively (e.g.
    `read -rs CONTROLLER_PASSWORD`) so it never appears on screen or in shell
    history.

# Part 1 — Build the stack

To build a developer "stack" from the fragalysis-backend repository we build the
fragalysis-backend container image first, and then switch to the fragalysis-stack
project where we build and push the final application image. We don't build the
front-end and instead rely on the frontend that is picked up in the stack Dockerfile
(e.g. xchem/fragalysis-frontend:latest).

`STACK_NAMESPACE` was already collected above; the remaining build variables are
derived and `export`ed before building: -

- BE_NAMESPACE is always "xchem"
- BE_IMAGE_TAG is the current fragalysis-backend branch name *slugified*
  (e.g. "m2ms-1234"). This is the `GITHUB_REF_SLUG` CI uses, so it is lower-cased
  and has `/` replaced by `-` (e.g. branch `feature/Foo` becomes `feature-foo`).
  Use the slug, not the raw branch name, so the tag matches any CI-built image.
- STACK_IMAGE_TAG is the same as BE_IMAGE_TAG
- FE_NAMESPACE is "xchem" (we don't build the front-end)
- FE_IMAGE_TAG is "latest" (the front-end image picked up by the stack Dockerfile)

With `STACK_NAMESPACE` already exported, set up the derived variables with: -

```
export BE_NAMESPACE=xchem
export BE_IMAGE_TAG=$(git rev-parse --abbrev-ref HEAD | tr '[:upper:]' '[:lower:]' | tr '/' '-')
export STACK_IMAGE_TAG=${BE_IMAGE_TAG}
export FE_NAMESPACE=xchem
export FE_IMAGE_TAG=latest
```

>   If the developer is on staging or production there's no need to build anything,
    as an official image will exist in `xchem/fragalysis-stack` DockerHub registry
    (typically `xchem/fragalysis-stack:latest`)

>   If there are no local modifications a build might not be necessary,
    as the CI process (`build-dev.yaml`) ensures that a container image is built for
    the current branch. If you're on branch `m2ms-1234` for example you might find
    the image `xchem/fragalysis-backend:m2ms-1234` already exists.

To build the fragalysis-backend AMD image: -

```
docker buildx build . --platform linux/amd64 --load \
  -t ${BE_NAMESPACE}/fragalysis-backend:${BE_IMAGE_TAG}
```

>   The `--load` is important. The stack Dockerfile does
    `FROM ${BE_NAMESPACE}/fragalysis-backend:${BE_IMAGE_TAG}`, so the image must be
    in the *local* image store. Without `--load` (or `--push`) a `docker-container`
    or `cloud` builder leaves the image only in its build cache, and the stack build
    silently pulls the published image from DockerHub instead of your local one -
    so your local backend changes would be lost. (`--load` works here because this
    is a single-platform build.)

Once the build has been successful we then move to the fragalysis-stack's project
directory.

>   If it is not present you can clone it from https://github.com/xchem/fragalysis-stack
    and use the `master` branch: -

```
cd ../fragalysis-stack
```

The stack image is pushed to the developer's own namespace, so log in first
(if you are not already): -

```
docker login
```

Then we can build and push the fragalysis-stack image: -

```
docker buildx build . --platform linux/amd64 \
  -t ${STACK_NAMESPACE}/fragalysis-stack:${STACK_IMAGE_TAG} \
  --build-arg BE_NAMESPACE=${BE_NAMESPACE} --build-arg BE_IMAGE_TAG=${BE_IMAGE_TAG} \
  --build-arg FE_NAMESPACE=${FE_NAMESPACE} --build-arg FE_IMAGE_TAG=${FE_IMAGE_TAG} \
  --push
```

>   `FE_NAMESPACE` and `FE_IMAGE_TAG` must be set (see the environment variables
    above). If they are unset the `--build-arg` expands to an empty value, which
    overrides the Dockerfile defaults (`xchem`/`latest`) with an empty string and
    produces an invalid `FROM /fragalysis-frontend:` reference.

# Part 2 — Deploy the stack

With the stack image built and pushed we can deploy it to Kubernetes. We trigger a
deployment through **AWX** (Ansible AWX / Automation Controller). Each developer has
their own AWX **Job Template** that, when launched, deploys their personal Fragalysis
Stack to the cluster. This part launches that Job Template and waits for it to finish.

The AWX credentials (`CONTROLLER_HOST`, `CONTROLLER_USERNAME`, `CONTROLLER_PASSWORD`)
were already gathered and exported up front, so the `awx` CLI is ready to use.

## 1. Derive the Job Template name

The Job Template is named `User (<User>) Developer Fragalysis Stack`, where
`<User>` is the developer's AWX username with its **first letter capitalised**
(the rest of the username is left as-is). For example `alan` becomes `Alan`, so
the template is `User (Alan) Developer Fragalysis Stack`.

Derive it from `CONTROLLER_USERNAME` so the two can never drift apart: -

```
USER_CAPITALISED="$(tr '[:lower:]' '[:upper:]' <<< "${CONTROLLER_USERNAME:0:1}")${CONTROLLER_USERNAME:1}"
JOB_TEMPLATE="User (${USER_CAPITALISED}) Developer Fragalysis Stack"
```

## 2. Launch the Job Template

We launch the template with the `awx` CLI (from `alanbchristie-awxkit`, an
unofficial fork of Ansible's `awxkit` that works on Python 3.13 and 3.14 — see the
install note below). The Job Template needs
two **extra variables** telling it which stack image to deploy — these identify the
image built and pushed in Part 1, so they are derived from the same variables:

- `stack_image` is `${STACK_NAMESPACE}/fragalysis-stack`
- `stack_image_tag` is `${STACK_IMAGE_TAG}`

`--monitor` streams the job output and makes the command exit non-zero if the
deployment fails, so an error is never swallowed: -

```
awx job_templates launch "${JOB_TEMPLATE}" \
  --extra_vars "stack_image: ${STACK_NAMESPACE}/fragalysis-stack
stack_image_tag: ${STACK_IMAGE_TAG}" \
  --monitor --wait
```

>   `--extra_vars` takes YAML (or JSON). The two-line value above is YAML; keep the
    keys unquoted and let the shell expand the variables. If you prefer JSON,
    `--extra_vars '{"stack_image": "...", "stack_image_tag": "..."}'` works too.

>   If the `awx` command is not installed, install it with
    `pip install 'alanbchristie-awxkit~=1.0'` — an unofficial fork of Ansible's
    `awxkit` maintained for Python 3.13 and 3.14. It installs the same `awxkit`
    package and `awx` command as the original, so don't install it alongside the
    upstream `awxkit` distribution (uninstall that first if it is present).

>   The launch resolves the template by name. If AWX reports that the name is
    ambiguous or not found, list the developer's templates to confirm the exact
    name with `awx job_templates list --all` and check the capitalisation of
    `<User>`.

A successful run reports the job as `successful`; the developer's Fragalysis Stack
is then deploying (or redeploying) on Kubernetes.
