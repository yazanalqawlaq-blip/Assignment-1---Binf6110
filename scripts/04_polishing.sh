#!/usr/bin/env bash
#SBATCH --job-name=polish
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=slurm/polish_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

ASSEMBLY="output_files/flye_assembly/assembly.fasta"
FASTQ="SRR32410565.fastq.gz"
MODEL="r1041_e82_400bps_sup_v5.0.0"
THREADS=16

[[ -f "$FASTQ" ]]    || { echo "FASTQ not found" >&2; exit 1; }
[[ -f "$ASSEMBLY" ]] || { echo "Assembly not found" >&2; exit 1; }

module load apptainer/1.2.4

apptainer exec medaka_latest.sif medaka_consensus -m "$MODEL" \
  -i "$FASTQ" \
  -d "$ASSEMBLY" \
  -o output_files/medaka_out \
  -t "$THREADS"

