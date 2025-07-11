#!/bin/bash

#SBATCH -p serial
#SBATCH --mem=8g
#SBATCH -J job_enzyme_adjustment_test_new_2000

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "enzyme_adjustment_test_new_2000;quit"

echo Job running on compute node `uname -n`
