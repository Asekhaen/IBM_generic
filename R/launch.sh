#!/bin/bash

############################################################
# LAUNCH
#
# Launch pipeline from command line.
#
# Command line usage:
#   sh launch.sh
#
############################################################

# Define R version
r_version="4.4.1"

# module load gcc/12.2.0  # NOTE: This is pre-loaded on Setonix
module load r/${r_version}

# Call main R launch script
Rscript launch.R

