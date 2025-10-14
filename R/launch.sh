#!/bin/bash

# You can edit the values of these directives as needed
#SBATCH --job-name=MyRJob         # Name of the job
#SBATCH --partition=work          # Partition to submit the job to
#SBATCH --ntasks=1                # Number of tasks (use 1 for a serial R job)
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --time=02:00:00           # Maximum wall-clock time (HH:MM:SS)
#SBATCH --mem=4G                  # Memory requested per node
#SBATCH --output=r_job_%j.out     # Standard output file
#SBATCH --error=r_job_%j.err      # Standard error file

# These are your executable commands, which you must leave unhashed.
# Load the desired R module, after confirming the version exists on Setonix.
module load r/4.4.1

# Navigate to the directory where your scripts are located
cd $SLURM_SUBMIT_DIR

# Execute the R script using Rscript
Rscript run_file.R
