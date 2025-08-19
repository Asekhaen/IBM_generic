#!/bin/bash

#SBATCH --job-name=vaccine3vs4
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

############################################################
# SUBMIT
#
# Submit array of cluster tasks.
#
# Written by A.J.Shattock
############################################################

# Extract inputs
job_type=$1
job_total=$2
n_cores=$3
hpc=$4
r_version=$5
log_file=$6

# ---- Load R ----

# Identify which cluster we're working from
case $hpc in

  # Load R on SciCore
  *scicore*)
    module purge
    module load R/${r_version}-foss-2023b
    ;;
    
  # Load R (and dependencies) on Pawsey-setonix
  *setonix*)
    module load gcc-native/12.3
    module load r/${r_version}
    ;;
esac

# ---- Prepare jobs to run ----

# Function to be called by srun
run_task() {
  
  # Calculate the base task number for this array job 
  job_base=$(( (SLURM_ARRAY_TASK_ID - 1) * n_cores ))
  
  # The actual task number is base + local task ID
  job_id=$(( job_base + SLURM_PROCID + 1 ))

  # Only run if task number is <= total jobs
  if (( job_id <= job_total )); then
  
    # Run R script which calls appropriate simulation function
    Rscript submit.R $job_type $job_id
    
    # Once complete, write job ID to log file
    echo $job_id >> $log_file
  fi
}

# Export function and all input arguments
export -f run_task
export n_cores job_total job_type log_file

# ---- Launch Jobs ----

# Call srun to run n_cores times in parallel
srun --ntasks=$n_cores bash -c "run_task"

