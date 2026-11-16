#!/usr/bin/env bash
# Submit from pipeline repo:
#   mkdir -p SLURM_LOGS
#   sbatch 20261112_launch_r84180_20261110.sbatch
# Or run the all-samples sheet:
#   sbatch --export=ALL,INPUT_FILE=20261112_samplesheet_r84180_20261110_all.yaml 20261112_launch_r84180_20261110.sbatch
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=batch
#SBATCH --output=SLURM_LOGS/.slurm_r84180_20261110_%j.out
#SBATCH --error=SLURM_LOGS/.slurm_r84180_20261110_%j.err
#SBATCH --mail-type=ALL
#SBATCH --job-name=ISA-r84180-20261110

set -euo pipefail
cd /Users/markpampuch/Downloads/tmp/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS

# module load nextflow
# module load singularity/3.9.7

export TMOUT=172800
export NXF_OPTS='-Xms3G -Xmx5G'
# export NXF_SINGULARITY_CACHEDIR=/ibex/scratch/projects/c2303/NXF_SINGULARITY_CACHEDIR
# export NXF_APPTAINER_CACHEDIR=/ibex/scratch/projects/c2303/NXF_APPTAINER_CACHEDIR
# export NXF_WORK=/ibex/scratch/projects/c2303/work

INPUT_FILE="/Users/markpampuch/Downloads/tmp/20260726_NEXTFLOW_INTEGRATION_SITE_ANALYSIS/20261116_samplesheet_r84180_20261110_LOCAL-TESTS.yaml"
# CONFIG_FILES="nextflow.config,conf/20261112_r84180_20261110.config"
CONFIG_FILES="nextflow.config"
NXF_OUTPUT_DIR="OUTPUTS/r84180_20261110_$(date -Iseconds | sed 's/-//g; s/://g; s/T/_/; s/+.*//')"
NXF_LOG_FILE="$NXF_OUTPUT_DIR/nextflow.log"

mkdir -p "$NXF_OUTPUT_DIR" SLURM_LOGS
cp "$INPUT_FILE" "$NXF_OUTPUT_DIR/"
# cp nextflow.config conf/20261112_r84180_20261110.config "$NXF_OUTPUT_DIR/"
cp nextflow.config "$NXF_OUTPUT_DIR/"
cp "$0" "$NXF_OUTPUT_DIR/slurm_script.sbatch"
echo "$(pwd)" > "$NXF_OUTPUT_DIR/projectDir.txt"

nextflow -log "$NXF_LOG_FILE" run . \
    -c nextflow.config \
    -with-report \
    -profile docker \
    --input "$INPUT_FILE" \
    --outdir "$NXF_OUTPUT_DIR" \
    --integration_dbs combined \
    -resume
