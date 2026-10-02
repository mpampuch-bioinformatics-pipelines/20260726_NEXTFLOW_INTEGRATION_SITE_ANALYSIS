# Respect `.cursorignore` Rules

Always adhere to the ignore patterns defined in `.cursorignore` across all tools (file viewing, directory listing, searching, background operations, and subagent delegation).

## Key Ignore Principles

1. **Large Data & Binary Files**:
   - Do NOT attempt to view or search large sequencing/genomic data files (`*.fastq`, `*.fastq.gz`, `*.fasta`, `*.fa`, `*.bam`, `*.sam`, `*.vcf`, `*.bed`, `*.gtf`, `*.gff`, `*.bw`, `*.h5`, `*.loom`, `*.mtx`, etc.) unless located under `TEST_DATA/`.
   - Do NOT view binary/archive files (`*.zip`, `*.tar`, `*.gz`, `*.sif`, `*.img`, `*.pdf`, `*.png`, `*.jpg`).
   - Do NOT inspect raw text/csv/tsv data dumps (`*.csv`, `*.tsv`, `*.txt`, `*.dat`, `*.raw`) EXCEPT for samplesheets matching `*samplesheet*` or test files in `TEST_DATA/`.

2. **Ignored Directories**:
   - Exclude Nextflow work and cache directories (`work/`, `WORK/`, `cache/`, `CACHE/`, `.nextflow/`, `.lineage/`, `work_symlink/`).
   - Exclude output, build, and temporary directories (`results/`, `output/`, `logs/`, `temp/`, `tmp/`, `scratch/`, `DATA/`, `data/`, `raw_data/`, `processed_data/`, `dist/`, `build/`).
   - Exclude environment and dependency folders (`venv/`, `env/`, `myenv/`, `.mamba*`, `.conda*`, `node_modules/`, `__pycache__/`).

3. **Explicitly Allowed Exceptions**:
   - **Samplesheets**: Files matching `*samplesheet*.csv`, `*samplesheet*.tsv`, `*samplesheet*.txt` are allowed.
   - **Logs**: Specific output/slurm log files in `OUTPUTS/**/*.log*`, `OUTPUTS/**/*.out*`, `OUTPUTS/**/*.err*`, `SLURM_LOGS/**/*.log*`, `SLURM_LOGS/**/*.out*`, `SLURM_LOGS/**/*.err*` are allowed.
   - **Test Data**: All files located inside `TEST_DATA/` are allowed.
