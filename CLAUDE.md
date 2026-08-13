# Claude Instructions for csgenetics_scrnaseq

## THIS REPOSITORY IS PUBLIC
Everything committed here is world-readable the moment it is pushed, and stays readable
in the git history even after deletion. Before committing, check that you are not adding:
- private object-storage paths (literal bucket names must be unambiguously public) or private/internal
  infrastructure identifiers; documented, reviewed, intentionally public code-controlled IDs are permitted
- customer/site/run identifiers (e.g. `EXAMPLE_SITE_RUN_001`) or colleague names
- credentials, credential file paths, or anything under `agent-secrets/`
- local machine paths (`/nssd2/...`)
- planning docs, handoffs, campaign logs, or other agent working files

Write agent working files to the session scratchpad or a private location, never into the
tree. Ad-hoc sample sheets that point at internal data belong in access-controlled storage
outside this repo. The `check-no-internal-files` CI job enforces the mechanical parts of
this. It scans the Git index and rejects non-public literal buckets in `s3`, `s3a`, `s3n`
URIs and standard AWS S3 HTTPS endpoints. In CI the index is the committed tip; locally it
is the content staged for the next commit. It also rejects Conda recipes whose local
source path resolves outside this repository. This is not an all-history scan, and
anything already published remains published and must be remediated or rotated as
appropriate.

## When developing
- Never write code that fails silently. Always fail loud.
- Never run the pipeline unless I explicitly ask you to do so.
- Never tell me 'I'm absolutely right' when I'm not right.
- Critically appraise what I say and don't assume that I'm right but do always do as I say.
- Don't do defensive programming. Trust properly implemented upstream solutions and don't add redundant safety checks downstream.

## When doing pull request reviews
- Always be sure to check the file content rather than just the diff.
- Before claiming duplicates or conflicts, provide explicit line numbers with actual values to verify the claim.
- When suggesting code simplifications, test all edge cases (null, false, undefined, empty) match the original behavior.
- Before recommending alternatives, verify they don't violate the same constraints you're trying to fix.

## When troubleshooting an issue
- Don't search the work directory for the work dir of a process by simply searching for files.
- Rather, always search the latest .nextflow.log file to look for the problematic process, find the work dir and go and search it that way.
- When troubleshooting an error, by default you should check the work dir of the failed process and check the logs to see what happened.

## When invalidating Nextflow cache for a specific process
- The .nextflow.log shows partial directory hashes (e.g., [eb/2ee0b1]). Always use wildcards to delete the full directory: rm -rf work/eb/2ee0b1*
- Search for BOTH "Cached" and "Submitted" entries: grep "process_name" .nextflow.log | grep -E "\[[a-f0-9]{2}/[a-f0-9]{6}\]"
- First show the user the full list of directories to be deleted, then delete them all in a single clear command.
