#!/usr/bin/env bash
#SBATCH --job-name=align_variants
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=slurm/align_variants_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

POLISHED="output_files/medaka_out/consensus.fasta"
REF="input_data/reference_fasta.fna"
READS="SRR32410565.fastq.gz"
OUTDIR="output_files"

[[ -f "$READS" ]]     || { echo "FASTQ not found" >&2; exit 1; }
[[ -f "$POLISHED" ]]  || { echo "Assembly not found" >&2; exit 1; }
[[ -f "$REF" ]]       || { echo "Reference not found" >&2; exit 1; }

module load StdEnv/2023 minimap2/2.28 samtools/1.20 bcftools/1.22

echo "Aligning reads to reference..."
minimap2 -a -x map-ont -t 16 "$REF" "$READS" | \
  samtools view -bS - | \
  samtools sort -@ 16 -o "$OUTDIR/reads_to_ref.sorted.bam" -

samtools index "$OUTDIR/reads_to_ref.sorted.bam"

# Variant calling with bcftools
# Using -Q 10 (min base quality) and -q 5 (min mapping quality) to reduce
# noise from low-confidence basecalls and poorly mapped reads
echo "Calling variants..."
bcftools mpileup -f "$REF" -Ou \
  -Q 10 -q 5 -d 8000 \
  --annotate FORMAT/DP,FORMAT/AD,INFO/AD \
  "$OUTDIR/reads_to_ref.sorted.bam" | \
bcftools call -m --ploidy 1 -Ov -o "$OUTDIR/raw_calls.vcf"

# Normalize variant representation before filtering
bcftools norm -f "$REF" "$OUTDIR/raw_calls.vcf" \
  -Ov -o "$OUTDIR/normalized_calls.vcf"

# Apply quality and depth filters
bcftools filter -e 'QUAL<30 || FORMAT/DP<15' \
  "$OUTDIR/normalized_calls.vcf" -Ov -o "$OUTDIR/filtered.vcf"

# Extract only variant sites (exclude reference-matching positions)
bcftools view -v snps,indels "$OUTDIR/filtered.vcf" \
  -Ov -o "$OUTDIR/variants_only.vcf"

bcftools stats "$OUTDIR/variants_only.vcf" > "$OUTDIR/variant_stats.txt"

echo "Aligning assembly to reference..."
minimap2 -a -x asm5 -t 16 "$REF" "$POLISHED" | \
  samtools view -bS - | \
  samtools sort -@ 16 -o "$OUTDIR/assembly_to_ref.sorted.bam" -

samtools index "$OUTDIR/assembly_to_ref.sorted.bam"

echo "Alignment and variant calling complete!"

# --- Summary ---
# Aligns raw reads to the reference with minimap2 (map-ont), calls
# variants with bcftools (mpileup/call), normalizes and filters them,
# then aligns the polished assembly to the reference for
# structural comparison.
#
# References:
# minimap2: https://github.com/lh3/minimap2
# samtools: https://www.htslib.org/doc/samtools.html
# bcftools: https://samtools.github.io/bcftools/bcftools.html
