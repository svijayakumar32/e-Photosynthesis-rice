#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=8g
#SBATCH -J job_array_gpmain_rice_420
#SBATCH -a 1-15

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "gpmain_rice_420_new;quit"

echo This is job task ${SLURM_ARRAY_TASK_ID} running on compute node `uname -n`