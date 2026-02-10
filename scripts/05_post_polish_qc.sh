#!/usr/bin/env bash
#SBATCH --job-name=post_polish_qc
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/post_polish_qc_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

POLISHED="output_files/medaka_out/consensus.fasta"
REF="input_data/reference_fasta.fna"

[[ -f "$POLISHED" ]] || { echo "Polished assembly not found" >&2; exit 1; }
[[ -f "$REF" ]]      || { echo "Reference not found" >&2; exit 1; }

module load apptainer/1.2.4

echo "Running QUAST on polished assembly..."
apptainer exec quast_5.2.0--py39pl5321h2add14b_1.sif quast.py \
  "$POLISHED" \
  -r "$REF" \
  -o output_files/quast_polished \
  --threads 8

echo "Running BUSCO on polished assembly..."
apptainer exec busco_5.7.1--pyhdfd78af_0.sif busco \
  -i "$POLISHED" \
  -o output_files/busco_polished \
  -m genome \
  --auto-lineage-prok \
  -c 8

echo "Post-polish QC complete!"

# --- Summary ---
# Repeats QUAST and BUSCO on the Medaka-polished assembly to
# quantify improvements in accuracy and completeness relative
# to the raw assembly.
#
# References:
# QUAST: https://github.com/ablab/quast
# BUSCO: https://busco.ezlab.org/
