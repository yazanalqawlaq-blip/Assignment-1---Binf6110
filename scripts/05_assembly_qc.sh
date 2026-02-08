#!/usr/bin/env bash
#SBATCH --job-name=assembly_qc
#SBATCH --account=def-cottenie
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/05_qc_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Assembly quality assessment ==="
echo "Started: $(date)"

ASSEMBLY_RAW="output_files/flye_assembly/assembly.fasta"
ASSEMBLY_POLISHED="output_files/medaka_polished/consensus.fasta"
REFERENCE="input_data/reference_fasta.fna"
THREADS=8

module load apptainer/1.2.4

# Pull containers
if [[ ! -f quast.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1
fi

if [[ ! -f busco.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/busco:5.7.1--pyhdfd78af_0
fi

if [[ ! -f bandage.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/bandage:0.8.1--hc9558a2_2
fi

# QUAST on raw assembly
echo "Running QUAST on raw assembly..."
apptainer exec quast.sif quast.py \
    ${ASSEMBLY_RAW} \
    -r ${REFERENCE} \
    -o output_files/quast_raw \
    --threads ${THREADS}

# QUAST on polished assembly
echo "Running QUAST on polished assembly..."
apptainer exec quast.sif quast.py \
    ${ASSEMBLY_POLISHED} \
    -r ${REFERENCE} \
    -o output_files/quast_polished \
    --threads ${THREADS}

# BUSCO on polished assembly
echo "Running BUSCO..."
apptainer exec busco.sif busco \
    -i ${ASSEMBLY_POLISHED} \
    -o output_files/busco_polished \
    -m genome \
    --auto-lineage-prok \
    -c ${THREADS}

# Bandage visualization
echo "Creating assembly graph visualization..."
apptainer exec bandage.sif Bandage image \
    output_files/flye_assembly/assembly_graph.gfa \
    figures/assembly_graph.png \
    --height 2000 --width 2000

# Copy key results
cp output_files/quast_polished/report.txt output_files/polished_quast_report.txt
cp output_files/quast_polished/report.pdf output_files/polished_quast_report.pdf
cp output_files/busco_polished/short_summary*.txt output_files/polished_busco_summary.txt

echo "=== QC complete: $(date) ==="
