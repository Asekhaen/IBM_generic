#!/bin/bash

#SBATCH --account=pawsey1163

#SBATCH --job-name=SDM.SA

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --time=23:00:00

#SBATCH --mem=200G

#SBATCH -o ./output%j.%a.out

#SBATCH --mail-user=andrew.coates@curtin.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_90

module load r/4.4.1

cd $MYSCRATCH/PSHB

Rscript src/spartan_call.R