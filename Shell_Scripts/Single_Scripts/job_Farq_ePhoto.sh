#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=5g
#SBATCH -J job_Farq_ePhoto

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "Farq_ePhoto_compariosn;quit"

echo Job running on compute node `uname -n`