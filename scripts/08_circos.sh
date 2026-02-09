#!/usr/bin/env bash
#SBATCH --job-name=circos_plot
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm/circos_plot_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

POLISHED="output_files/medaka_out/consensus.fasta"
REFERENCE="input_data/reference_fasta.fna"

[[ -f "$POLISHED" ]] || { echo "Assembly not found" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found" >&2; exit 1; }

module load StdEnv/2023 scipy-stack/2024b

echo "Creating genome comparison visualization..."

python3 << 'PYEOF'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ref_data = []
with open("input_data/reference_fasta.fna") as f:
    seq = ""
    for line in f:
        if line.startswith(">"):
            if seq:
                ref_data.append(len(seq))
                seq = ""
        else:
            seq += line.strip()
    if seq:
        ref_data.append(len(seq))

asm_data = []
with open("output_files/medaka_out/consensus.fasta") as f:
    seq = ""
    for line in f:
        if line.startswith(">"):
            if seq:
                asm_data.append(len(seq))
                seq = ""
        else:
            seq += line.strip()
    if seq:
        asm_data.append(len(seq))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

colors_ref = ['#3498db', '#2ecc71', '#9b59b6', '#f39c12']
colors_asm = ['#e74c3c', '#1abc9c', '#34495e']

y_pos = 0
for i, length in enumerate(ref_data):
    ax1.barh(y_pos, length, color=colors_ref[i % len(colors_ref)], 
             edgecolor='black', linewidth=0.5, label=f'Sequence {i+1}: {length:,} bp')
    y_pos += 1

ax1.set_yticks(range(len(ref_data)))
ax1.set_yticklabels([f'Seq {i+1}' for i in range(len(ref_data))])
ax1.set_xlabel('Length (bp)', fontsize=12, fontweight='bold')
ax1.set_title('Reference Genome (S. enterica LT2)', fontsize=14, fontweight='bold')
ax1.legend(loc='best', fontsize=9)
ax1.grid(axis='x', alpha=0.3)

y_pos = 0
for i, length in enumerate(asm_data):
    ax2.barh(y_pos, length, color=colors_asm[i % len(colors_asm)], 
             edgecolor='black', linewidth=0.5, label=f'Contig {i+1}: {length:,} bp')
    y_pos += 1

ax2.set_yticks(range(len(asm_data)))
ax2.set_yticklabels([f'Contig {i+1}' for i in range(len(asm_data))])
ax2.set_xlabel('Length (bp)', fontsize=12, fontweight='bold')
ax2.set_title('Assembled Genome (Polished)', fontsize=14, fontweight='bold')
ax2.legend(loc='best', fontsize=9)
ax2.grid(axis='x', alpha=0.3)

plt.suptitle('Genome Structure Comparison', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('output_files/genome_comparison.png', dpi=300, bbox_inches='tight')
print("SUCCESS: Genome comparison plot saved!")
PYEOF

echo "Visualization complete!"
ls -lh output_files/genome_comparison.png
