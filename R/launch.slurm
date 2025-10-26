#!/bin/bash

# You can edit the values of these directives as needed
#SBATCH --account=pawsey1163
#SBATCH --job-name=parameter_check  # Name of the job
#SBATCH --nodes=1
#SBATCH --ntasks=1                # Number of tasks (use 1 for a serial R job)
#SBATCH --cpus-per-task=4         # Number of CPUs per task
#SBATCH --time=20:00:00           # Maximum wall-clock time (HH:MM:SS)
#SBATCH --mem-per-cpu=4G          # Memory requested per node
#SBATCH --output=r_job_%j.out     # Standard output file
#SBATCH --error=r_job_%j.err      # Standard error file

#SBATCH --mail-user=john.ohiolei@postgrad.curtin.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_90

# These are your executable commands, which you must leave unhashed.
# Load the desired R module, after confirming the version exists on Setonix...
module load r/4.4.1

# Navigate to the directory where your scripts are located
cd $MYSOFTWARE/IBM_generic/R

# Execute the R script using Rscript
Rscript run_file.R
