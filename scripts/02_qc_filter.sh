#!/usr/bin/env bash
#SBATCH --job-name=qc_filter
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=slurm/02_qc_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Quality control and filtering ==="
echo "Started: $(date)"

ACC="SRR32410565"
READS="${ACC}.fastq.gz"
THREADS=8

module load apptainer/1.2.4

# Pull containers
if [[ ! -f nanoplot.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_0
fi

if [[ ! -f nanofilt.sif ]]; then
    apptainer pull docker://quay.io/biocontainers/nanofilt:2.8.0--py_0
fi

# Run NanoPlot QC
echo "Running NanoPlot..."
apptainer exec nanoplot.sif NanoPlot \
    --fastq ${READS} \
    --outdir output_files/nanoplot_raw \
    --threads ${THREADS}

# Filter reads: Q≥10, length≥1000bp
echo "Filtering reads (Q≥10, length≥1000bp)..."
apptainer exec nanofilt.sif gunzip -c ${READS} | \
    NanoFilt -q 10 -l 1000 | \
    gzip > ${ACC}.filtered.fastq.gz

# QC on filtered reads
echo "QC on filtered reads..."
apptainer exec nanoplot.sif NanoPlot \
    --fastq ${ACC}.filtered.fastq.gz \
    --outdir output_files/nanoplot_filtered \
    --threads ${THREADS}

echo "=== QC complete: $(date) ==="
