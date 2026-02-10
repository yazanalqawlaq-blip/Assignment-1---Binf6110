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

RAW_ASM="output_files/flye_assembly/assembly.fasta"
READS="SRR32410565.fastq.gz"
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"
THREADS=16

[[ -f "$READS" ]]   || { echo "FASTQ not found" >&2; exit 1; }
[[ -f "$RAW_ASM" ]] || { echo "Assembly not found" >&2; exit 1; }

module load apptainer/1.2.4

apptainer exec medaka_latest.sif medaka_consensus -m "$MEDAKA_MODEL" \
  -i "$READS" \
  -d "$RAW_ASM" \
  -o output_files/medaka_out \
  -t "$THREADS"

# --- Summary ---
# Polishes the raw assembly with Medaka using the nanopore reads
# and the R10.4.1 super-accuracy basecalling model to correct
# remaining consensus errors.
#
# References:
# Medaka: https://github.com/nanoporetech/medaka
