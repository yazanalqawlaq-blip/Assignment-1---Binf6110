#!/usr/bin/env bash
#SBATCH --job-name=alignment
#SBATCH --account=def-cottenie
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=slurm/06_alignment_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Alignment to reference ==="
echo "Started: $(date)"

ACC="SRR32410565"
READS="${ACC}.filtered.fastq.gz"
ASSEMBLY="output_files/medaka_polished/consensus.fasta"
REFERENCE="input_data/reference_fasta.fna"
THREADS=16

module load StdEnv/2023 minimap2/2.28 samtools/1.19

# Alignment 1: Raw reads to reference (for variant calling)
echo "Aligning reads to reference..."
minimap2 -ax map-ont \
    -t ${THREADS} \
    ${REFERENCE} \
    ${READS} | \
    samtools view -bS - | \
    samtools sort -@ ${THREADS} -o output_files/reads_to_ref.sorted.bam -

samtools index output_files/reads_to_ref.sorted.bam

# Alignment 2: Assembly to reference (for structural comparison)
echo "Aligning assembly to reference..."
minimap2 -ax asm5 \
    -t ${THREADS} \
    ${REFERENCE} \
    ${ASSEMBLY} | \
    samtools view -bS - | \
    samtools sort -@ ${THREADS} -o output_files/assembly_to_ref.sorted.bam -

samtools index output_files/assembly_to_ref.sorted.bam

# Calculate alignment statistics
echo "Calculating alignment statistics..."
samtools flagstat output_files/reads_to_ref.sorted.bam > output_files/reads_alignment_stats.txt
samtools flagstat output_files/assembly_to_ref.sorted.bam > output_files/assembly_alignment_stats.txt

echo "=== Alignment complete: $(date) ==="
