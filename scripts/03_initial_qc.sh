#!/usr/bin/env bash
#SBATCH --job-name=initial_qc
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/initial_qc_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

ASSEMBLY="output_files/flye_assembly/assembly.fasta"
REFERENCE="input_data/reference_fasta.fna"

[[ -f "$ASSEMBLY" ]]  || { echo "Assembly not found" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found" >&2; exit 1; }

module load apptainer/1.2.4

echo "Running QUAST on raw assembly..."
apptainer exec quast_5.2.0--py39pl5321h2add14b_1.sif quast.py \
  "$ASSEMBLY" \
  -r "$REFERENCE" \
  -o output_files/quast_raw \
  --threads 8

echo "Running BUSCO on raw assembly..."
apptainer exec busco_5.7.1--pyhdfd78af_0.sif busco \
  -i "$ASSEMBLY" \
  -o output_files/busco_raw \
  -m genome \
  --auto-lineage-prok \
  -c 8

echo "Initial QC complete!"
