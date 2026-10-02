# Use nf-core modules before writing a process

Do not write a Nextflow process from scratch when nf-core already publishes that tool.

Before adding a module:

1. `module load tower-cli`
2. `nf-core modules list remote` and find the tool.
3. If it exists, install it with `nf-core modules add MODULE/PATH` (for example `seqkit/stats`, `busco/download`).
4. Write a local module only when the remote list has no module that does the job.

Pipeline-specific logic that no remote module covers (clustering BLAST hits, extracting a flanked coordinate range) stays local. Prefer a remote module plus a thin local step over reimplementing the tool.
