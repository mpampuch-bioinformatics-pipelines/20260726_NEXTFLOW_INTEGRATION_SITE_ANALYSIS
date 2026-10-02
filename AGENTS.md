# Workspace Guidelines and Rules

## `.cursorignore` Compliance
All operations, file reads, searches, and subagent tasks must strictly respect the `.cursorignore` rules:

- **Ignore Large & Raw Data**: Do not read or search large data files (`*.fastq`, `*.fasta`, `*.bam`, `*.vcf`, `*.bed`, `*.bw`, `*.h5`, `*.loom`, `*.mtx`), generic `*.csv`/`*.tsv`/`*.txt`/`*.dat`/`*.raw` files, or binary/archive/image assets (`*.pdf`, `*.png`, `*.jpg`, `*.zip`, `*.tar`, `*.sif`).
- **Ignore Work & Output Directories**: Skip Nextflow runtime dirs (`work/`, `.nextflow/`, `.lineage/`, `cache/`), data dirs (`data/`, `DATA/`), and output/build/venv dirs (`results/`, `output/`, `tmp/`, `venv/`, `node_modules/`, `__pycache__/`).
- **Allowed Exceptions**:
  - Samplesheets (`*samplesheet*.csv`, `*samplesheet*.tsv`, `*samplesheet*.txt`)
  - Log files under `OUTPUTS/` and `SLURM_LOGS/` (`*.log`, `*.out`, `*.err`)
  - Any test files located inside `TEST_DATA/`

Refer to [.cursorignore](file:///ibex/scratch/projects/c2303/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/.cursorignore) and [.agents/rules/cursorignore.md](file:///ibex/scratch/projects/c2303/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/.agents/rules/cursorignore.md) for complete details.
