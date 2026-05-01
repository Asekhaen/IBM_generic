#!/bin/bash

#SBATCH --account=pawsey1163
#SBATCH --job-name=allee_effect_load   # Name of the job
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1                # Number of tasks (use 1 for a serial R job)
#SBATCH --cpus-per-task=4        # Number of CPUs per task
#SBATCH --time=5:00:00           # Maximum wall-clock time (HH:MM:SS)
#SBATCH --mem-per-cpu=4G          # Memory requested per node
#SBATCH --output=r_job_%j.out     # Standard output file
#SBATCH --error=r_job_%j.err      # Standard error file

#SBATCH --mail-user=john.ohiolei@postgrad.curtin.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_90

# These are your executable commands, which you must leave unhashed.
# Load the desired R module, after confirming the version exists on Setonix.
module load r/4.4.1

# Navigate to the directory where your scripts are located
cd $MYSOFTWARE/Allee_effect/R

# Execute the R script using Rscript
Rscript run_file.R