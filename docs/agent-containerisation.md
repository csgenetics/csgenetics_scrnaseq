# Autonomous Agent Containerisation

> **Public repository.** This document is visible to the world. It contains no
> secrets, no host-specific paths, and no IAM identifiers. All deployment-specific
> values are injected via environment variables at runtime.
>
> This containerised-agent pattern was developed by CS Genetics for its internal
> repositories and applied here as the third rollout under epic GRT-352. The
> pattern is the same across repos; this document is self-contained so that
> external readers can understand the agent's architecture and security model
> without access to CS Genetics' private repositories.

This document describes how the `csgenetics-scrnaseq-agent` Claude Code agent
operates against this repository: its architecture, security model, setup,
runtime behaviour, and rotation procedures.

The goal is to let a Claude Code agent work on this repository autonomously with
a **bounded blast radius**: the agent can create branches, commit, push, and open
PRs (targeting `devel`), but it cannot merge to `main` or `devel`, cannot touch
any other GitHub repo, cannot access secrets or files on the host outside
tightly-scoped bind mounts, and cannot escape the container.

---

## Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [Security Model](#security-model)
- [Setup](#setup)
- [Running the Agent](#running-the-agent)
- [Authentication](#authentication)
- [Rotation Procedures](#rotation-procedures)
- [Troubleshooting](#troubleshooting)

---

## Overview

Each invocation of the agent:

1. Starts a fresh Docker container from the `scrnaseq-agent:latest` image
2. Bind-mounts this repository at `/scrnaseq_git_repo` (read-write)
3. Bind-mounts scoped AWS credentials at `/home/node/.aws` (read-only)
4. Bind-mounts the GitHub App private key at `/run/secrets/github-app.pem` (read-only)
5. Injects the Claude Max subscription OAuth token as an environment variable
6. Runs `claude -p` with Opus 4.6 in 1M-context mode and `--dangerously-skip-permissions`
7. Exits and auto-removes the container (via `--rm`)

Agent-authored PRs target **`devel`**, not `main`. `devel` is the integration
branch; `main` receives release promotions from `devel`. Both branches are
protected — the bot can merge neither.

The relevant files, all under `docker/agent/`:

| File | Runs on | Purpose |
|---|---|---|
| `Dockerfile` | build-time | Defines the image (Node 22, git, gh, AWS CLI, Claude Code) |
| `entrypoint.sh` | container | Configures git identity and credential helper, runs Claude |
| `gh-credential-helper.sh` | container | Mints fresh GitHub App installation tokens on demand |
| `gh-shim.sh` | container | Wraps `gh` so every invocation gets a fresh token |
| `run-on-beast.sh` | host | Validates secrets, bind-mounts, launches the container |

---

## Architecture

### Five conceptual components

**1. Your host.** Holds the long-lived secrets and the repo clone. The directory
layout is determined by `$REPO_PATH` and `$SECRETS_DIR` (environment variables
you set before running the wrapper script). A typical layout:

```
$SECRETS_DIR/
├── claude-oauth-token              # Claude Max subscription OAuth token
├── github/
│   └── scrnaseq-agent-app.pem      # GitHub App private key
└── aws/
    ├── credentials                 # Scoped IAM access key
    └── config                      # Region and output format
```

All files should be mode 600; the secrets directory itself mode 700.

**2. The container image (`scrnaseq-agent:latest`).** Built from
`docker/agent/Dockerfile`. Runs as uid 1000 (the `node` user from the
`node:22-slim` base image). Contains `git`, `gh`, `jq`, `openssl`, `curl`,
`unzip`, AWS CLI v2, and `@anthropic-ai/claude-code`.

**3. The host-side wrapper (`docker/agent/run-on-beast.sh`).** Validates that all
required files exist on the host, loads the Claude OAuth token, and launches
`docker run` with the right mounts and env vars.

**4. The GitHub App (`csgenetics-scrnaseq-agent`).** Org-owned, installed only on
`csgenetics/csgenetics_scrnaseq`, with `contents:write`, `pull_requests:write`,
`issues:write`, and nothing else. The app's private key lives on the host and is
mounted read-only into the container.

**5. Branch protection on `main` AND `devel`.** Both branches are configured so
only allowlisted users can push or merge. The bot is not on the allowlist. This
is the primary defence that prevents the bot from merging its own PRs — enforced
server-side by GitHub, regardless of what permissions the installation token
claims to have.

### Data flow

```
┌────────────────────────┐
│ operator invokes        │
│ run-on-beast.sh "..."   │
└────────────┬────────────┘
             │
             ▼
┌──────────────────────────────────┐
│ host-side wrapper                │
│ - reads claude-oauth-token       │
│ - validates all secrets exist    │
│ - exec docker run with mounts   │
└────────────┬─────────────────────┘
             │
             ▼
┌──────────────────────────────────────────────┐
│ container                                    │
│  entrypoint.sh → git config → claude -p ...  │
│  git push → credential helper → JWT → token  │
│  gh pr create → shim → fresh token → gh      │
└──────────────────────────────────────────────┘
             │
             ▼
       container exits, JSON result on stdout
```

---

## Security Model

### What the container CAN do

- Read and write files in `/scrnaseq_git_repo`
- Read scoped AWS credentials from `/home/node/.aws`
- Mint and use short-lived GitHub installation tokens (1 hour)
- Make outbound network requests
- Create feature branches, commit, push, open PRs (targeting `devel`)
- Access S3 buckets the scoped IAM user is authorised for

### What the container CANNOT do, by design

| Attempted action | Defence |
|---|---|
| Merge a PR into `main` or `devel` | Branch protection `restrictions` allowlist |
| Push directly to `main` or `devel` | Branch protection |
| See any other GitHub repo | GitHub App installation scoping |
| Read host files outside the bind mounts | Container filesystem namespacing |
| Read the host's personal AWS credentials | Only scoped agent credentials are mounted |
| Escape the container | Standard Docker isolation |
| Persist anything outside the mounted repo | Container removed on exit (`--rm`) |
| Bypass branch protection as admin | `enforce_admins: true` |

### Why `--dangerously-skip-permissions` is safe here

The container itself is the isolation boundary. Claude's per-tool permission
prompts are bypassed, but the container limits what files and network resources
are accessible. A "rogue Claude" inside the container can at worst write junk to
the mounted repo, which is cleaned up with `git checkout .`.

---

## Setup

### Prerequisites

- A Linux host with Docker 27+ and the invoking user in the `docker` group
- `git`, `gh`, `openssl`, `jq`, `curl` on the host
- The invoking user's uid should match the in-container `node` user (uid 1000)
- A clone of `csgenetics/csgenetics_scrnaseq`
- A Claude Max subscription (for the OAuth token)
- Admin access to the csgenetics GitHub org

### Step 1: Create a scoped AWS IAM user

Create an IAM user (e.g. `scrnaseq-agent`) with an inline policy scoped to only
the S3 buckets the agent needs. Template policy:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "ReadOnlyBucketLevelOps",
      "Effect": "Allow",
      "Action": ["s3:ListBucket", "s3:GetBucketLocation", "s3:ListBucketMultipartUploads"],
      "Resource": [
        "arn:aws:s3:::csgx.cbk",
        "arn:aws:s3:::csgx.graticule",
        "arn:aws:s3:::themis.bios.grt",
        "arn:aws:s3:::csg-nextflow",
        "arn:aws:s3:::csg-reference",
        "arn:aws:s3:::csg-tower-bucket",
        "arn:aws:s3:::csgx-circleci",
        "arn:aws:s3:::csgx.external.data",
        "arn:aws:s3:::csgx.latch",
        "arn:aws:s3:::csgx.public.readonly",
        "arn:aws:s3:::long-term-sequencing-data"
      ]
    },
    {
      "Sid": "ReadOnlyObjects",
      "Effect": "Allow",
      "Action": ["s3:GetObject", "s3:GetObjectVersion"],
      "Resource": [
        "arn:aws:s3:::csgx.cbk/*",
        "arn:aws:s3:::csgx.graticule/*",
        "arn:aws:s3:::themis.bios.grt/*",
        "arn:aws:s3:::csg-nextflow/*",
        "arn:aws:s3:::csg-reference/*",
        "arn:aws:s3:::csg-tower-bucket/*",
        "arn:aws:s3:::csgx-circleci/*",
        "arn:aws:s3:::csgx.external.data/*",
        "arn:aws:s3:::csgx.latch/*",
        "arn:aws:s3:::csgx.public.readonly/*",
        "arn:aws:s3:::long-term-sequencing-data/*"
      ]
    },
    {
      "Sid": "WriteNfTestFixtures",
      "Effect": "Allow",
      "Action": [
        "s3:PutObject", "s3:DeleteObject",
        "s3:AbortMultipartUpload", "s3:ListMultipartUploadParts"
      ],
      "Resource": "arn:aws:s3:::csg-reference/internal_nf_tests_data/*"
    }
  ]
}
```

Generate an access key for the user and store the credentials on the host:

```
$SECRETS_DIR/aws/credentials    # [default] profile with access key
$SECRETS_DIR/aws/config          # [default] profile with region
```

Both files mode 600.

### Step 2: Obtain a Claude Max subscription OAuth token

```bash
claude setup-token
```

Copy the printed token into `$SECRETS_DIR/claude-oauth-token` (mode 600). Use an
editor, not shell redirection, to keep the token out of shell history.

### Step 3: Create the GitHub App

In the csgenetics org settings, create a GitHub App:

- **Name:** `csgenetics-scrnaseq-agent`
- **Webhook:** unchecked
- **Repository permissions:** Contents (R/W), Pull requests (R/W), Issues (R/W)
- **All other permissions:** No access
- **Install:** only on `csgenetics/csgenetics_scrnaseq`

Note the **App ID** and **Installation ID**. Generate a private key and store it
at `$SECRETS_DIR/github/scrnaseq-agent-app.pem` (mode 600).

### Step 4: Configure branch protection

Apply identical protection rules to **both `main` AND `devel`**:

```bash
for BRANCH in main devel; do
  gh api --method PUT \
    -H "Accept: application/vnd.github+json" \
    /repos/csgenetics/csgenetics_scrnaseq/branches/$BRANCH/protection \
    --input - <<'JSON'
{
  "required_status_checks": null,
  "enforce_admins": true,
  "required_pull_request_reviews": {
    "dismiss_stale_reviews": false,
    "require_code_owner_reviews": false,
    "required_approving_review_count": 0,
    "require_last_push_approval": false
  },
  "restrictions": {
    "users": ["didillysquat", "emilyscher"],
    "teams": [],
    "apps": []
  },
  "required_linear_history": false,
  "allow_force_pushes": false,
  "allow_deletions": false,
  "required_conversation_resolution": true,
  "lock_branch": false,
  "allow_fork_syncing": false
}
JSON
done
```

Verify both branches:

```bash
for B in main devel; do
  echo "=== $B ===";
  gh api /repos/csgenetics/csgenetics_scrnaseq/branches/$B/protection | jq '{
    allowlist_users: [.restrictions.users[].login],
    enforce_admins: .enforce_admins.enabled
  }';
done
```

### Step 5: Build the image

```bash
cd $REPO_PATH
docker build -t scrnaseq-agent:latest -f docker/agent/Dockerfile .
```

### Step 6: First run

```bash
REPO_PATH=/path/to/csgenetics_scrnaseq \
SECRETS_DIR=/path/to/agent-secrets/scrnaseq \
./docker/agent/run-on-beast.sh
```

With no prompt argument, the entrypoint runs a diagnostic that reports the
working directory, git identity, and first few directory entries as JSON.

---

## Running the Agent

```bash
REPO_PATH=/path/to/csgenetics_scrnaseq \
SECRETS_DIR=/path/to/agent-secrets/scrnaseq \
./docker/agent/run-on-beast.sh "Your prompt here"
```

PRs should always target **`devel`**, not `main`.

The entrypoint hardcodes `--model claude-opus-4-6[1m]` and `--output-format json`.
Override the model by passing a flag after the prompt:

```bash
./docker/agent/run-on-beast.sh "quick question" --model claude-sonnet-4-6
```

### Environment variable overrides

| Env var | Default | Purpose |
|---|---|---|
| `REPO_PATH` | *(required)* | Host path to the repo checkout |
| `SECRETS_DIR` | *(required)* | Parent dir with secrets |
| `DOCKER_IMAGE` | `scrnaseq-agent:latest` | Image tag |
| `GITHUB_APP_ID` | `3409188` | GitHub App ID (public identifier) |
| `GITHUB_APP_INSTALLATION_ID` | `124688983` | Installation ID (public identifier) |

---

## Authentication

Three independent authentication systems are active inside the container:

**Claude Max subscription.** OAuth token injected as `CLAUDE_CODE_OAUTH_TOKEN` env
var. Never stored on disk inside the container.

**GitHub App installation tokens.** The RSA private key is mounted read-only at
`/run/secrets/github-app.pem`. The credential helper signs a JWT and mints a
fresh 1-hour installation access token on every git/gh invocation. Tokens are
never stored. Scope: only `csgenetics/csgenetics_scrnaseq`.

**AWS scoped IAM user.** Credentials at `/home/node/.aws` (bind-mounted
read-only). Scope is determined by the IAM policy attached to the user — see
[Setup Step 1](#step-1-create-a-scoped-aws-iam-user) for the template.

---

## Rotation Procedures

### Claude OAuth token (yearly)

Regenerate with `claude setup-token`, replace `$SECRETS_DIR/claude-oauth-token`.
No container restart needed — `run-on-beast.sh` reads the file fresh on every run.

### GitHub App private key (every ~6 months)

Generate a new key in the GitHub App settings. Replace
`$SECRETS_DIR/github/scrnaseq-agent-app.pem`. Delete the old key in the GitHub
UI after verifying the new one works.

### AWS credentials (yearly or on personnel change)

Create a new access key in IAM. Replace `$SECRETS_DIR/aws/credentials`.
Deactivate and later delete the old key.

---

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `run-on-beast.sh` errors before container starts | Missing file in `$SECRETS_DIR` — the wrapper validates all inputs |
| Claude reports "Not logged in" | `CLAUDE_CODE_OAUTH_TOKEN` env var missing or token expired |
| `git push` fails with 401/403 | Installation token minting failed — check `openssl rsa -check` on the `.pem` file |
| `git push` to `main` or `devel` returns 403 | Expected — branch protection working correctly |
| `gh pr merge` returns 405 | Expected — bot is not on the branch protection allowlist |
| Permission errors writing to `/scrnaseq_git_repo` | Container uid (1000) doesn't match host file ownership |
| `aws: command not found` | Image needs rebuild (`docker build ...`) — AWS CLI is installed in the Dockerfile |

---

*This document was created in April 2026 as part of the containerised-agent rollout (GRT-355, epic GRT-352). Keep it up to date when you change any of the referenced files or procedures.*
