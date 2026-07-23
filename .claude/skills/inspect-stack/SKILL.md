---
name: inspect-stack
description: >-
  Download and inspect a deployed Fragalysis Stack's backend log
  (logs/backend.log) from its Pod running in Kubernetes. The KUBECONFIG is
  supplied by name when the skill is run or via the KUBECONFIG environment
  variable, and the user names the Namespace the stack belongs to. Use when
  the user says "inspect the stack", "get the stack log", "check the backend
  log", or "/inspect-stack <kubeconfig> <namespace>".
---

# Inspect a Fragalysis Stack

This skill retrieves the backend log from a **deployed** Fragalysis Stack so it
can be inspected. The stack runs in **Kubernetes**, so we reach it with
`kubectl`, locate the stack **Pod** in the user's **Namespace**, copy its
`logs/backend.log` out of the container, and then summarise what the log shows.

We need two things before we can talk to the cluster: a **KUBECONFIG** (which
cluster) and a **Namespace** (which stack). Gather both up front.

# 1. Resolve the KUBECONFIG

The KUBECONFIG can arrive two ways — accept whichever the user provides:

- **By name when the skill is run** — the user passes a path, e.g.
  `/inspect-stack nw-xch-dev.yaml my-namespace`. A named file takes precedence.
- **Via the `KUBECONFIG` environment variable** — a KUBECONFIG file already
  named in the environment.

Resolve to a single file, confirm it exists, and export it so every `kubectl`
call in this skill uses it. Do **not** print the file's contents — it holds a
cluster credential.

```bash
# KUBECONFIG_ARG is the path the user named when running the skill, if any.
if [ -n "${KUBECONFIG_ARG}" ]; then
  export KUBECONFIG="${KUBECONFIG_ARG}"
elif [ -n "${KUBECONFIG}" ]; then
  export KUBECONFIG="${KUBECONFIG}"
else
  echo "No KUBECONFIG: name one when running the skill or set the KUBECONFIG env var." >&2
  exit 1
fi

if [ ! -f "${KUBECONFIG}" ]; then
  echo "KUBECONFIG file '${KUBECONFIG}' does not exist." >&2
  exit 1
fi
echo "Using KUBECONFIG=${KUBECONFIG}"
```

>   The Fragalysis Rancher clusters' KUBECONFIG files are produced by the
    `get-rancher-kubeconfigs` skill (e.g. `nw-xch-dev.yaml`, `nw-xch-prod.yaml`).
    If the user has none, point them there first.

# 2. Confirm the Namespace

The **Namespace** identifies which stack to inspect and is **required**. It is
named when the skill is run (the second argument above) — if it was not
provided, ask the user for it rather than guessing. Verify it exists on the
cluster so a typo fails clearly instead of later commands returning nothing:

```bash
# NAMESPACE is the namespace the user named.
if [ -z "${NAMESPACE}" ]; then
  echo "No Namespace given: the user must name the Namespace the stack belongs to." >&2
  exit 1
fi
if ! kubectl get namespace "${NAMESPACE}" >/dev/null 2>&1; then
  echo "Namespace '${NAMESPACE}' not found on this cluster (check the KUBECONFIG and the name)." >&2
  exit 1
fi
```

# 3. Find the stack Pod

The backend log lives inside the **stack Pod** in that Namespace. Find the
running Pod:

```bash
POD=$(kubectl -n "${NAMESPACE}" get pods \
  --field-selector=status.phase=Running \
  -o jsonpath='{.items[*].metadata.name}')
echo "Pod(s): ${POD}"
```

Usually the Namespace holds a single stack Pod, so `POD` is that one name. If
more than one name comes back, list the Pods so the user can see them and pick
the stack Pod (the one running the backend); do not guess when it is ambiguous:

```bash
kubectl -n "${NAMESPACE}" get pods
```

If no Pod is returned the stack is not running (or is still starting) — report
that and stop.

# 4. Download `logs/backend.log`

Copy the backend log out of the stack Pod to a local file so we can inspect it
without holding the container open. Read it through a shell (`sh -c`) so the
container's working directory applies and the relative `logs/backend.log` path
resolves. Save the download into the scratchpad directory:

```bash
LOCAL_LOG="${TMPDIR:-/tmp}/backend-${NAMESPACE}.log"
kubectl -n "${NAMESPACE}" exec "${POD}" -- sh -c 'cat logs/backend.log' > "${LOCAL_LOG}"
echo "Downloaded backend log to ${LOCAL_LOG}"
wc -l "${LOCAL_LOG}"
```

>   If the Pod runs more than one container, `kubectl` will ask you to name one —
    re-run the `exec` with `-c <container>` for the stack (backend) container.
    You can list a Pod's containers with
    `kubectl -n "${NAMESPACE}" get pod "${POD}" -o jsonpath='{.spec.containers[*].name}'`.

>   If `logs/backend.log` is not found at the relative path, the container's
    working directory is not the code root. The repo is served from `/code`, so
    fall back to the absolute path:
    `kubectl -n "${NAMESPACE}" exec "${POD}" -- cat /code/logs/backend.log`.

# 5. Inspect and report

Now inspect the downloaded log and give the user a useful summary rather than
dumping the whole file. Prefer the dedicated Read/Grep tools on `${LOCAL_LOG}`:

- Show the **tail** — the most recent activity is usually what matters.
- Surface **errors and warnings** (e.g. `ERROR`, `CRITICAL`, `Traceback`,
  `Exception`, `WARNING`) with enough surrounding context to be actionable.
- Note the **time range** the log covers and anything that looks like a
  startup, migration, or Celery/task problem.

Report the local path, a short summary of what the log shows, and quote the
specific lines that back up any problem you call out.
