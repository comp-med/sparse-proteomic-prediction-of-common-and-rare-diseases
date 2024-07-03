#!/usr/bin/bash

#SBATCH --job-name=dz_prot_pred_expansion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=100:00:00
#SBATCH --array=4,6,8,11,15,24,28,39,43,45,46,49,51,169,187
#SBATCH --mail-user=jc050021@gsk.com
#SBATCH --mail-type=ALL
#SBATCH --output=logs/%x_%j_report.out
#SBATCH --error=logs/%x_%j_log.err

# load module
module add R/4.0.2

# asign directories
work_dir="/home/ukbb_inc_dz_prediction/Explore_1536_Expansion_analyses/02_inc_dz_prediction/bin"
log_dir="$work_dir/logs"
script="$work_dir/01_ML_1536_pipeline_6m_protein_pxinfo_pred_subc.R"

cd ${work_dir}
# input parameter file 
export f_input="../data_input/${1}"

#run as job array
echo "Job ID: $SLURM_ARRAY_TASK_ID"
dz="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var{print $1}' ${f_input})"
date="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var{print $2}' ${f_input})"
yrs="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var{print $3}' ${f_input})"

echo ${dz} ${date} ${yrs} 

R CMD BATCH --no-save --no-restore "--args ${dz} ${date} ${yrs}" $script $log_dir/res_1536_${dz}_${yrs}.Rout

echo "Done."

