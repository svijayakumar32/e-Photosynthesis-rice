#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=8g
#SBATCH -J job_array_gpmain_rice_all
#SBATCH -a 1-10

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "gpmain_rice_new(140);quit"

echo This is job task ${SLURM_ARRAY_TASK_ID}
