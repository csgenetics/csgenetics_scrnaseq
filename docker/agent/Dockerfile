# syntax=docker/dockerfile:1.6
#
# Container image for the csgenetics-scrnaseq-agent autonomous bot.
#
# Contents:
#   - Node 22 (for Claude Code CLI)
#   - git, gh CLI (for repo operations)
#   - jq, curl, openssl, unzip (for scripting and JWT signing)
#   - AWS CLI v2 (for S3 access — credentials bind-mounted at /home/node/.aws)
#   - @anthropic-ai/claude-code (the agent runtime)
#
# Runs as non-root user `node` (uid 1000). If your host user has a different
# uid, bind-mount permissions may need adjustment.
#
# Secrets are NEVER baked into the image.  At `docker run` time the caller:
#   - Mounts the GitHub App private key read-only at /run/secrets/github-app.pem
#   - Mounts scoped AWS credentials read-only at /home/node/.aws
#   - Injects CLAUDE_CODE_OAUTH_TOKEN, GITHUB_APP_ID, GITHUB_APP_INSTALLATION_ID
#     as env vars
#
# Build from the repo root:
#   docker build -t scrnaseq-agent:latest -f docker/agent/Dockerfile .

FROM node:22-slim

# Base system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates \
      curl \
      git \
      gnupg \
      jq \
      openssl \
      unzip \
    && rm -rf /var/lib/apt/lists/*

# gh CLI (official Debian package repo)
RUN install -d -m 0755 /etc/apt/keyrings \
 && curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
      | tee /etc/apt/keyrings/githubcli-archive-keyring.gpg > /dev/null \
 && chmod go+r /etc/apt/keyrings/githubcli-archive-keyring.gpg \
 && echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" \
      > /etc/apt/sources.list.d/github-cli.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends gh \
 && rm -rf /var/lib/apt/lists/*

# AWS CLI v2 (official installer).
# Architecture-aware: node:22-slim is multi-arch, so resolve the AWS CLI
# archive name from dpkg's arch rather than hardcoding x86_64.
RUN case "$(dpkg --print-architecture)" in \
      amd64) AWSCLI_ARCH=x86_64 ;; \
      arm64) AWSCLI_ARCH=aarch64 ;; \
      *) echo "Unsupported architecture: $(dpkg --print-architecture)" >&2; exit 1 ;; \
    esac \
 && curl -sSL "https://awscli.amazonaws.com/awscli-exe-linux-${AWSCLI_ARCH}.zip" -o /tmp/awscliv2.zip \
 && unzip -q /tmp/awscliv2.zip -d /tmp \
 && /tmp/aws/install \
 && rm -rf /tmp/awscliv2.zip /tmp/aws \
 && aws --version

# Claude Code CLI
# TODO: pin to a specific version for reproducibility
RUN npm install -g @anthropic-ai/claude-code && claude --version

# Shorthand for `claude --dangerously-skip-permissions`. Safe inside this
# container because the container is the blast-radius boundary (see
# docs/agent-containerisation.md).
RUN printf '#!/bin/sh\nexec claude --dangerously-skip-permissions "$@"\n' \
      > /usr/local/bin/claude-yolo \
 && chmod 0755 /usr/local/bin/claude-yolo

# Entrypoint and helpers (owned by root, world-executable, read-only for agent)
COPY docker/agent/entrypoint.sh           /usr/local/bin/entrypoint.sh
COPY docker/agent/gh-credential-helper.sh /usr/local/bin/gh-credential-helper.sh
COPY docker/agent/gh-shim.sh              /usr/local/bin/gh
RUN chmod 0755 /usr/local/bin/entrypoint.sh \
               /usr/local/bin/gh-credential-helper.sh \
               /usr/local/bin/gh

# Non-root user.  The node:22-slim base image provides a 'node' user at
# uid 1000.  If your host user's uid differs, bind-mount permissions will
# break — adjust accordingly.
USER node
WORKDIR /scrnaseq_git_repo

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
