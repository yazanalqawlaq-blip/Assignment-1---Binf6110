#!/usr/bin/env bash
#SBATCH --job-name=polish
#SBATCH --account=def-cottenie
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/04_polish_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Medaka polishing ==="
echo "Started: $(date)"

ACC="SRR32410565"
ASSEMBLY="output_files/flye_assembly/assembly.fasta"
READS="${ACC}.filtered.fastq.gz"
OUTDIR="output_files/medaka_polished"
MODEL="r1041_e82_400bps_sup_v500"
THREADS=16

module load apptainer/1.2.4

# Pull Medaka container
if [[ ! -f medaka.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/medaka:2.0.1--py311h87f3376_0
fi

# Run Medaka consensus
echo "Polishing with Medaka model: ${MODEL}..."
apptainer exec medaka.sif medaka_consensus \
    -i ${READS} \
    -d ${ASSEMBLY} \
    -o ${OUTDIR} \
    -m ${MODEL} \
    -t ${THREADS}

echo "Polished assembly: ${OUTDIR}/consensus.fasta"
echo "=== Polishing complete: $(date) ==="
