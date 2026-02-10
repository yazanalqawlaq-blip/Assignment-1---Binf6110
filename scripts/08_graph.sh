#!/usr/bin/env bash
#SBATCH --job-name=graph
#SBATCH --account=def-cottenie
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=slurm/graph_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

cd /scratch/yazanalq/Assignment-1---Binf6110

module load StdEnv/2023 scipy-stack/2024b

python3 << 'EOF'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

os.makedirs("output_files", exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 6))
metrics = ['N50\n(Mb)', 'Total Length\n(Mb)', 'BUSCO\nCompleteness (%)', 'Mismatches\nper 100 kbp', 'Indels\nper 100 kbp']
raw = [3.32, 5.11, 99.2, 27.41, 3.83]
polished = [3.32, 5.11, 98.4, 27.11, 3.71]

x = np.arange(len(metrics))
width = 0.35

bars1 = ax.bar(x - width/2, raw, width, label='Raw Assembly', edgecolor='black', linewidth=1)
bars2 = ax.bar(x + width/2, polished, width, label='Polished Assembly', edgecolor='black', linewidth=1)

ax.set_ylabel('Value', fontsize=12, fontweight='bold')
ax.set_title('Assembly Quality Metrics: Raw vs Polished', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=10)
ax.legend(fontsize=10)
ax.grid(axis='y', alpha=0.3, linestyle='--')

for bars in (bars1, bars2):
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h, f'{h:.2f}',
                ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('output_files/figure3_quality_graph.png', dpi=300, bbox_inches='tight')
print("Figure saved: output_files/figure3_quality_graph.png")
EOF

ls -lh output_files/figure3_quality_graph.png

# --- Summary ---
# Produces a grouped bar chart comparing key quality metrics
# (N50, total length, BUSCO, mismatches, indels) between the
# raw and polished assemblies using matplotlib.
#
# References:
# matplotlib: https://matplotlib.org/stable/index.html
