#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=5g
#SBATCH -J job_gpmain_rice_340

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "gpmain_rice_340_new;quit"

echo Job running on compute node `uname -n`