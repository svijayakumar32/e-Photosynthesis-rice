#!/bin/bash

#SBATCH -p test
#SBATCH -J testjob_twotimes

source /etc/profile
module add matlab/2022a

matlab -nodisplay -nodesktop -r "twotimes('5');quit"

echo Job running on compute node `uname -n`