#!/usr/bin/env bash
#SBATCH --job-name=visualization
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/visualization_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

POLISHED="output_files/medaka_out/consensus.fasta"
REF="input_data/reference_fasta.fna"
ALN_BAM="output_files/assembly_to_ref.sorted.bam"

[[ -f "$POLISHED" ]] || { echo "Assembly not found" >&2; exit 1; }
[[ -f "$REF" ]]      || { echo "Reference not found" >&2; exit 1; }
[[ -f "$ALN_BAM" ]]  || { echo "Alignments not found" >&2; exit 1; }

module load StdEnv/2023 samtools/1.20

echo "Extracting alignment statistics..."
samtools depth "$ALN_BAM" > output_files/assembly_coverage.txt
samtools stats "$ALN_BAM" > output_files/assembly_alignment_stats.txt

echo "Creating assembly comparison summary..."
cat > output_files/assembly_summary.txt << EOF
Assembly vs Reference Comparison
Generated: $(date)

Raw Assembly Stats:
$(grep -E "Total length|N50|contigs" output_files/quast_raw/report.txt | head -5)

Polished Assembly Stats:
$(grep -E "Total length|N50|contigs" output_files/quast_polished/report.txt | head -5)

Alignment Coverage:
Average depth: $(awk '{sum+=$3; count++} END {if(count>0) print sum/count; else print 0}' output_files/assembly_coverage.txt)
Total aligned bases: $(samtools view -c "$ALN_BAM")
EOF

echo "Visualization data prepared!"

# --- Summary ---
# Extracts coverage depth and alignment statistics from the
# assembly-to-reference BAM, then generates a text summary
# comparing raw vs polished QUAST metrics.
#
# References:
# samtools: https://www.htslib.org/doc/samtools.html
