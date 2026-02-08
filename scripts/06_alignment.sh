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

MEDAKA_ASSEMBLY="output_files/medaka_out/consensus.fasta"
REFERENCE="input_data/reference_fasta.fna"
FASTQ="SRR32410565.fastq.gz"
OUTDIR="output_files"

[[ -f "$FASTQ" ]] || { echo "FASTQ not found" >&2; exit 1; }
[[ -f "$MEDAKA_ASSEMBLY" ]] || { echo "Assembly not found" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found" >&2; exit 1; }

module load StdEnv/2023 minimap2/2.28 samtools/1.20 bcftools/1.22

echo "Aligning reads to reference..."
minimap2 -a -x map-ont -t 16 "$REFERENCE" "$FASTQ" | \
  samtools view -bS - | \
  samtools sort -@ 16 -o "$OUTDIR/reads_to_ref.sorted.bam" -

samtools index "$OUTDIR/reads_to_ref.sorted.bam"

echo "Calling variants..."
bcftools mpileup -f "$REFERENCE" -Ou -Q 7 -q 0 -d 10000 \
  -a FORMAT/DP,FORMAT/AD "$OUTDIR/reads_to_ref.sorted.bam" | \
bcftools call -m --ploidy 1 -Ov -o "$OUTDIR/raw_calls.vcf"

bcftools filter -e 'QUAL<20 || FORMAT/DP<10' \
  "$OUTDIR/raw_calls.vcf" -Ov -o "$OUTDIR/filtered.vcf"

bcftools view -v snps,indels "$OUTDIR/filtered.vcf" \
  -Ov -o "$OUTDIR/variants_only.vcf"

bcftools stats "$OUTDIR/variants_only.vcf" > "$OUTDIR/variant_stats.txt"

echo "Aligning assembly to reference..."
minimap2 -a -x asm5 -t 16 "$REFERENCE" "$MEDAKA_ASSEMBLY" | \
  samtools view -bS - | \
  samtools sort -@ 16 -o "$OUTDIR/assembly_to_ref.sorted.bam" -

samtools index "$OUTDIR/assembly_to_ref.sorted.bam"

echo "Alignment and variant calling complete!"
