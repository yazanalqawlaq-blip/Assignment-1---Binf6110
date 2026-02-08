#!/usr/bin/env bash
#SBATCH --job-name=variant_calling
#SBATCH --account=def-cottenie
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/07_variants_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Variant calling with Clair3 ==="
echo "Started: $(date)"

BAM="output_files/reads_to_ref.sorted.bam"
REFERENCE="input_data/reference_fasta.fna"
OUTDIR="output_files/clair3_variants"
MODEL="r1041_e82_400bps_sup_v500"
THREADS=16

module load apptainer/1.2.4

# Pull Clair3 container
if [[ ! -f clair3.sif ]]; then
    apptainer pull docker://hkubal/clair3:v1.0.10
fi

# Run Clair3
echo "Running Clair3 variant calling..."
apptainer exec clair3.sif /opt/bin/run_clair3.sh \
    --bam_fn=${BAM} \
    --ref_fn=${REFERENCE} \
    --threads=${THREADS} \
    --platform=ont \
    --model_path=/opt/models/${MODEL} \
    --output=${OUTDIR} \
    --include_all_ctgs

# Filter and extract high-quality variants
module load StdEnv/2023 bcftools/1.19

echo "Filtering variants (QUAL≥20, DP≥10)..."
bcftools view \
    -i 'QUAL>=20 && DP>=10' \
    ${OUTDIR}/merge_output.vcf.gz | \
    bcftools view -v snps,indels > output_files/variants_filtered.vcf

# Extract variant statistics
echo "Extracting variant statistics..."
bcftools stats output_files/variants_filtered.vcf > output_files/variant_stats.txt

# Count variants by type
echo "Variant counts:" > output_files/variant_summary.txt
echo "SNPs: $(bcftools view -v snps output_files/variants_filtered.vcf | grep -v '^#' | wc -l)" >> output_files/variant_summary.txt
echo "Indels: $(bcftools view -v indels output_files/variants_filtered.vcf | grep -v '^#' | wc -l)" >> output_files/variant_summary.txt

echo "=== Variant calling complete: $(date) ==="
