#!/usr/bin/env bash
#SBATCH --job-name=assembly
#SBATCH --account=def-cottenie
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/03_assembly_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Flye assembly ==="
echo "Started: $(date)"

ACC="SRR32410565"
READS="${ACC}.filtered.fastq.gz"
OUTDIR="output_files/flye_assembly"
THREADS=16

module load apptainer/1.2.4

# Pull Flye container
if [[ ! -f flye.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0
fi

# Run Flye
echo "Running Flye with nano-hq preset..."
apptainer exec flye.sif flye \
    --nano-hq ${READS} \
    --genome-size 5m \
    --out-dir ${OUTDIR} \
    --threads ${THREADS}

echo "Assembly: ${OUTDIR}/assembly.fasta"
echo "=== Assembly complete: $(date) ==="
