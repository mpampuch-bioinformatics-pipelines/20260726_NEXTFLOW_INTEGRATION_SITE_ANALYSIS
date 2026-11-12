# Project Profile — Integration Site Analysis

## Purpose

- **Pipeline name:** nf-core/integrationsiteanalysis (`NFCORE_INTEGRATIONSITEANALYSIS`)
- **Scientific question / intended use:** Genomic integration-site analysis from long reads (construct mapping, locus identification, QC and reporting)
- **Supported input types:** Samplesheet (YAML/CSV) with long-read FASTQ paths; reference assembly; BLAST query constructs; plasmid GenBank
- **Expected outputs:** Assembly/QC summaries, BLAST hit tables, locus intervals, IGV reports, MultiQC report
- **Maintainers:** Mark Pampuch

## Scientific contract

- **Organism(s):** Project-dependent (reference supplied per run)
- **Reference genome/build:** User-provided FASTA (`--reference`)
- **Annotation release:** Lifted from reference where applicable
- **Sample metadata schema:** nf-core samplesheet + `nextflow_schema.json` / `assets/schema_input.json`
- **Important assumptions / exclusions:** Long-read inputs; plasmid and construct sequences provided for mapping context

## Execution

- **Nextflow version:** As on IBEX (`module load nextflow` via `env.sh`)
- **DSL:** DSL2
- **Container / environment strategy:** Singularity/Apptainer on IBEX (`-profile singularity`); Wave-built images where configured
- **Local test command:** `nextflow run . -profile test,singularity --outdir <OUTDIR>`
- **nf-test command:** `nf-test test` (tagged runs per `agent-resources/AGENTS.md`); full suite when ready
- **Slurm profile:** `launch-job.sbatch` / `launch-test.sbatch`
- **Approved partitions/account/QOS:** IBEX defaults in sbatch launchers (confirm before production)
- **Scratch and output locations:** `/ibex/scratch/projects/c2303/…`; `work_symlink` → shared work tree
- **Is job submission allowed without approval?** Default: **NO** (see `agent-resources/AGENTS.md`)

## Definition of done

See `agent-resources/instructions/nextflow.md`. Update this profile when scientific behavior or supported inputs change. After `./agent-resources/uninstall.sh`, replace this line with a project-local note.
