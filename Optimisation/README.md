This folder contains scripts associated with the optimisation of the e-Photosynthesis model.

The evolutionary algorithm is run using the following scripts:
- `average.m` calculates the average Vmax in each generation
- `mutate.m` generates variation in Vmax populations by modifying mutatePercentage
- `rankPop.m` ranks populations based on their CO2 uptake values
- `resizePop.m` resizes the population after each generation
- `stdev.m` calculates the standard deviation of Vmax in each generation
- `twoPoint.m` is a variation of `resizePop.m`, but it is not used in the current simulations

These files are called by scripts in the `gpmain` folder, which optimise the model at a specific CO2 concentration (130-420 umol mol^-1 at intervals of 10 umol mol^-1).
As the optimisation process is computationally intensive, it is recommended that users run the scripts using a high performance computing (HPC) cluster, allotting sufficient memory of 8GB per job.

For the purposes of this work, all optimisations were run on the Lancaster University High End Computing (HEC) cluster using basic Linux commands.
- For further documentation on the Lancaster University HEC, please refer to https://lancaster-hec.readthedocs.io/en/latest/index.html. 
- For general documentation about SLURM, please see: https://slurm.schedmd.com/.

The MATLAB scripts in gpmain (`gpmain_rice_xxx_new.m`) require corresponding job scripts (`job_gpmain_rice_xxx.sh`.sh) that are submitted from the login node to the SLURM scheduling system.
In our analysis, the batch job scripts are run as arrays - enabling the same program to be run multiple times.

Each job script is a text file containing the commands to be run along with the compute resources required to run them. 
* Prior to running a script, it may be necessary to convert line endings of the text files from DOS (Windows) to Unix, e.g.

`dos2unix job_array_gpmain_rice_xxx.sh`

Batch jobs are run on the HEC by creating a batch job script and submitting it to the system using the command `sbatch`. 
For example, a job script named `job_array_gpmain_rice_xxx.sh` can be submitted with the command:

`sbatch job_array_gpmain_rice_xxx.sh`
