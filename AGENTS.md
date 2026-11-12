# Agent instructions

Canonical rules live in [`agent-resources/AGENTS.md`](agent-resources/AGENTS.md). Read that file before substantive work, then the playbook it names. Paths inside it are relative to `agent-resources/`.

Completion reports use the structure in that file: Changed, Validation, Scientific/operational notes, Not verified, Reproduction.

## Project-specific rules

Additional rules in [`.agents/rules/`](.agents/rules/):

- [nf-core-modules.md](.agents/rules/nf-core-modules.md) — prefer `nf-core modules add` over bespoke processes when a match exists
- [cursorignore.md](.agents/rules/cursorignore.md) — data and directory boundaries for reads and searches (see also [`.cursorignore`](.cursorignore))

When this pipeline no longer needs agent scaffolding, remove the stub and bundle with `./agent-resources/uninstall.sh` (project rules in `.agents/` are kept).
