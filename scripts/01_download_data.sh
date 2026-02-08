#!/usr/bin/env bash
#SBATCH --job-name=download_data
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=slurm/01_download_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Downloading Salmonella enterica data ==="
echo "Started: $(date)"

# Variables
ACC="SRR32410565"
WORK_DIR="${HOME}/Assignment-1---Binf6110"
THREADS=8

cd ${WORK_DIR}

# Download raw reads
echo "Downloading SRA data..."
wget -O ${ACC}.sra "https://sra-pub-run-odp.s3.amazonaws.com/sra/${ACC}/${ACC}"

# Load apptainer module
module load apptainer/1.2.4

# Pull SRA-tools container if not present
if [[ ! -f sra-tools.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
fi

# Convert to FASTQ
echo "Converting to FASTQ..."
apptainer exec sra-tools.sif fasterq-dump \
    ${ACC}.sra \
    --threads ${THREADS} \
    --outdir .

# Compress
gzip -f ${ACC}.fastq

# Download reference genome
echo "Downloading reference genome..."
cd input_data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
gunzip GCF_000006945.2_ASM694v2_genomic.fna.gz
mv GCF_000006945.2_ASM694v2_genomic.fna reference_fasta.fna

echo "=== Download complete: $(date) ==="
