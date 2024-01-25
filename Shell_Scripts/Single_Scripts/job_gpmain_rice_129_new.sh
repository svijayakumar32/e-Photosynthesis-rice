#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=8g
#SBATCH -J job_gpmain_rice_129

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "gpmain_rice_129_new(129);quit"

echo Job running on compute node `uname -n`