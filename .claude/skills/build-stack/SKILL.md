---
name: build-stack
description: Build a developer stack, which consists of backend and stack container images
---

To build a developer "stack" from the fragalysis-backend repository we build the
fragalysis-backend container image first, and then switch to the fragalysis-stack
project where we build and push the final application image. We don't build the
front-end and instead rely on the frontend that is picked up in the stack Dockerfile
(e.g. xchem/fragalysis-frontend:latest).

Environment variables are used to control the build: -

- BE_NAMESPACE is always "xchem"
- BE_IMAGE_TAG is set to the current fragalysis-backend branch name (e.g. "m2ms-1234", the GITHUB_REF_SLUG)
- STACK_NAMESPACE is the developer's GitHub username (e.g. "alanbchristie")
- STACK_IMAGE_TAG is the same as BE_IMAGE_TAG

>   If the developer is on staging or production there's no need to build anything,
    as an official image will exist in `xchem/fragalysis-stack` DockerHub registry
    (typically `xchem/fragalysis-stack:latest`)

>   If there are no local modifications a built might not be necessary,
    as the CI process (`build-dev.yaml`) ensures that a container image is built for
    the current branch. If you're on branch `m2ms-1234` for example you might find
    the image `xchem/fragalysis-backend:m2ms-1234` already exists.

To build the fragalysis-backend AMD image: -

```
docker buildx build . --platform linux/amd64 \
  -t ${BE_NAMESPACE}/fragalysis-backend:${BE_IMAGE_TAG}
```

Once the build has been successful we then move to the fragalysis-stack's project
directory.

>   If it is not present you can clone it from https://github.com/xchem/fragalysis-stack
    and use the `master` branch: -

```
cd ../fragalysis-stack
```

Then we can build and push the fragalysis-stack image: -

```
docker buildx build . --platform linux/amd64 \
  -t ${STACK_NAMESPACE}/fragalysis-stack:${STACK_IMAGE_TAG} \
  --build-arg BE_NAMESPACE=${BE_NAMESPACE} --build-arg BE_IMAGE_TAG=${BE_IMAGE_TAG} \
  --build-arg FE_NAMESPACE=${FE_NAMESPACE} --build-arg FE_IMAGE_TAG=${FE_IMAGE_TAG} \
  --push
```
