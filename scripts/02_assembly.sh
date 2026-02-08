#!/usr/bin/env bash
#SBATCH --job-name=assembly
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/assembly_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

FASTQ="SRR32410565.fastq.gz"

[[ -f "$FASTQ" ]] || { echo "FASTQ not found"; exit 1; }

module load apptainer/1.2.4

apptainer exec flye_2.9.6--py311h2de2dd3_0.sif flye \
    --nano-hq "$FASTQ" \
    --out-dir output_files/flye_assembly \
    --threads 16
