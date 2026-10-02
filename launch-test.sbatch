#!/usr/bin/env bash
# Submit: mkdir -p SLURM_LOGS && sbatch launch-test.sbatch
#SBATCH --time=04:00:00                     # BUSCO still runs even though hifiasm is skipped
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4                   # Match conf/test.config local executor (4 cpus / 8.GB)
#SBATCH --ntasks=1
#SBATCH --mem=16G                           # Nextflow JVM (~5G) plus local process limit (8G)
#SBATCH --partition=batch
#SBATCH --output=SLURM_LOGS/.slurm_%j.out
#SBATCH --error=SLURM_LOGS/.slurm_%j.err
#SBATCH --mail-type=ALL
#SBATCH --job-name=TEST-ISA-ASSEMBLIES

# Load modules
module load nextflow
module load singularity/3.9.7

# Configure environment variables
export NXF_OPTS='-Xms3G -Xmx5G'
export NXF_SINGULARITY_CACHEDIR=/ibex/scratch/projects/c2303/NXF_SINGULARITY_CACHEDIR
export NXF_APPTAINER_CACHEDIR=/ibex/scratch/projects/c2303/NXF_APPTAINER_CACHEDIR
export NXF_WORK=/ibex/scratch/projects/c2303/work

INPUT_FILE="samplesheet.test.20260910.pre-assembled.yaml"
CONFIG_FILE="nextflow.config"
NXF_OUTPUT_DIR="TESTS/TEST_OUTPUTS/$(date -Iseconds | sed 's/-//g; s/://g; s/T/_/; s/+.*//')"
NXF_LOG_FILE="$NXF_OUTPUT_DIR/nextflow.log"

mkdir -p "$NXF_OUTPUT_DIR" SLURM_LOGS
cp "$INPUT_FILE" "$NXF_OUTPUT_DIR"
cp "$CONFIG_FILE" "$NXF_OUTPUT_DIR"
cp "$0" "$NXF_OUTPUT_DIR/slurm_script.sbatch"
echo "$(pwd)" > "$NXF_OUTPUT_DIR/projectDir.txt"

# test profile sets --input to samplesheet.test.yaml (user-provided primary/alternate, skip hifiasm)
nextflow -log "$NXF_LOG_FILE" run . \
    -c "$CONFIG_FILE" \
    -with-report \
    -profile kaust \
    --input "$INPUT_FILE" \
    --outdir "$NXF_OUTPUT_DIR" \
    -resume
