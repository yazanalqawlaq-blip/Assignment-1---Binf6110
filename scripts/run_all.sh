#!/usr/bin/env bash
# Master script to run entire pipeline
# Submit each job with dependency on previous job

set -euo pipefail

echo "Submitting Salmonella assembly pipeline..."

# Submit jobs with dependencies
JOB1=$(sbatch --parsable scripts/01_download_data.sh)
echo "Job 1 (Download): ${JOB1}"

JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} scripts/02_qc_filter.sh)
echo "Job 2 (QC/Filter): ${JOB2}"

JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} scripts/03_assembly.sh)
echo "Job 3 (Assembly): ${JOB3}"

JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} scripts/04_polish.sh)
echo "Job 4 (Polish): ${JOB4}"

JOB5=$(sbatch --parsable --dependency=afterok:${JOB4} scripts/05_assembly_qc.sh)
echo "Job 5 (Assembly QC): ${JOB5}"

JOB6=$(sbatch --parsable --dependency=afterok:${JOB4} scripts/06_alignment.sh)
echo "Job 6 (Alignment): ${JOB6}"

JOB7=$(sbatch --parsable --dependency=afterok:${JOB6} scripts/07_variant_calling.sh)
echo "Job 7 (Variants): ${JOB7}"

echo "All jobs submitted! Monitor with: squeue -u $USER"
