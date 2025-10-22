#!/bin/bash

#SBATCH --account=pawsey1163
#SBATCH --job-name=parameter_sampling  # Name of the job
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1                # Number of tasks (use 1 for a serial R job)
#SBATCH --cpus-per-task=32        # Number of CPUs per task
#SBATCH --time=23:00:00           # Maximum wall-clock time (HH:MM:SS)
#SBATCH --array=1-500            # run 1000 parameter sets in parallel
#SBATCH --mem-per-cpu=4G          # Memory requested per node
#SBATCH --output=r_job_%j.out     # Standard output file
#SBATCH --error=r_job_%j.err      # Standard error file

#SBATCH --mail-user=john.ohiolei@postgrad.curtin.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_90

# Available R module, this version should exists on Setonix.
module load r/4.4.1

# directory hosting my scripts 
cd $MYSOFTWARE/IBM_generic/R

# This execute the R script using Rscript
Rscript run_parallel.R



