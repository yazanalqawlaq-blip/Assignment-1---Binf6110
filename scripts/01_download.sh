#!/usr/bin/env bash
#SBATCH --job-name=download
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm/download_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

ACC="SRR32410565"

module load apptainer/1.2.4

wget -O ${ACC}.sra "https://sra-pub-run-odp.s3.amazonaws.com/sra/${ACC}/${ACC}"

apptainer exec sra-tools_3.2.1--h4304569_1.sif fasterq-dump ${ACC}.sra --threads 4

gzip ${ACC}.fastq

mkdir -p input_data
cd input_data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
gunzip GCF_000006945.2_ASM694v2_genomic.fna.gz
mv GCF_000006945.2_ASM694v2_genomic.fna reference_fasta.fna
