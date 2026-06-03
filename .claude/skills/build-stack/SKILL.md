---
name: build-stack
description: Build a developer stack, which consists of backend and stack container images
---

To build a developer "stack" from the fragalysis-backend repository we build the
fragalysis-backend container image first, and then switch to the fragalysis-stack
project where we build and push the final application image. We don't build the
front-end and instead rely on the frontend that is picked up in the stack Dockerfile
(e.g. xchem/fragalysis-frontend:latest).

Environment variables are used to control the build. **Ask the developer for only
one value, `STACK_NAMESPACE`** (their DockerHub/GitHub username, e.g.
"alanbchristie"); derive all the others as described below and `export` them before
building: -

- STACK_NAMESPACE is the developer's username — **ask for this**
- BE_NAMESPACE is always "xchem"
- BE_IMAGE_TAG is the current fragalysis-backend branch name *slugified*
  (e.g. "m2ms-1234"). This is the `GITHUB_REF_SLUG` CI uses, so it is lower-cased
  and has `/` replaced by `-` (e.g. branch `feature/Foo` becomes `feature-foo`).
  Use the slug, not the raw branch name, so the tag matches any CI-built image.
- STACK_IMAGE_TAG is the same as BE_IMAGE_TAG
- FE_NAMESPACE is "xchem" (we don't build the front-end)
- FE_IMAGE_TAG is "latest" (the front-end image picked up by the stack Dockerfile)

So, once you have `STACK_NAMESPACE` from the developer, set everything up with: -

```
export STACK_NAMESPACE=<the value the developer gave you>

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
