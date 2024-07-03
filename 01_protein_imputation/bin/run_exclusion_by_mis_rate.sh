#!/usr/bin/bash

#SBATCH --job-name=ukbb_exclusions
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=05:00:00
#SBATCH --mail-user=jc050021@gsk.com
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j_%u_report.out
#SBATCH --error=logs/%x_%j_%u_log.err

# load module
module add R/4.0.2

work_dir="/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/00_protein_imputation/bin"
log_dir="$work_dir/logs"
script="$work_dir/01_save_exlusion_by_missing_rate.R"

R CMD BATCH --no-save --no-restore $script $log_dir/Exclusions_by_mis_rate.Rout

echo "Done."

