#!/bin/bash

#SBATCH --job-name=report
#SBATCH --time=00:05:00
#SBATCH --mem=100MB
#SBATCH --cpus-per-task=1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

############################################################
# REPORT
#
# Submit array of cluster tasks.
#
# Written by A.J.Shattock
############################################################

# Extract inputs
array_id=$1
report_file=$2

# Write sacct job report to file
sacct -j $array_id --format=JobID,AllocCPUS,State,Elapsed -X > $report_file

