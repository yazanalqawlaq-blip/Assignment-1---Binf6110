#!/usr/bin/env bash
#SBATCH --job-name=visualizations
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/08_visualizations_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Creating visualizations ==="
echo "Started: $(date)"

module load StdEnv/2023 r/4.4.0

# Install required R packages if not present
Rscript -e "
packages <- c('ggplot2', 'dplyr', 'tidyr', 'vcfR', 'circlize', 'RColorBrewer', 'patchwork')
new_packages <- packages[!(packages %in% installed.packages()[,'Package'])]
if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')
"

# Run visualization script
Rscript scripts/08_visualizations.R

echo "=== Visualizations complete: $(date) ==="
