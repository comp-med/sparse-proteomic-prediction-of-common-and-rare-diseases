#!/usr/bin/bash

#SBATCH --job-name=ukbb_prot_impute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --exclusive
#SBATCH --time=192:00:00
#SBATCH --array=1-4
#SBATCH --mail-user=jc050021@gsk.com
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j_%u_report.out
#SBATCH --error=logs/%x_%j_%u_log.err

# load module
module add R/4.0.2

work_dir="/home/ukbb_inc_dz_prediction/03_inc_dz_prediction_1536_fix/bin"
log_dir="$work_dir/logs"
script="$work_dir/00_impute_protein_npx_Expansion.R"

cd ${work_dir}
# input parameter file 
export f_input="${1}"

#run as job array
echo "Job ID: $SLURM_ARRAY_TASK_ID"
panel="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var{print $1}' ${f_input})"

echo ${panel} 

R CMD BATCH --no-save --no-restore "--args ${panel}" $script $log_dir/update_impute_npx_missForest_${panel}.Rout

echo "Done."

