############################################################
# LAUNCH
#
# Launch pipeline for mosquito IBM model.
#
############################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

message("Running mosquito IBM")

# Load options into global environment
o = set_options(do_step = 2)  # See options.R

# 1) Prepare stuff
# run_prepare()  # See prepare.R

# 2) Run simulations
run_simulations()  # See simulations.R

# 3) Perform stats analysis
# run_analysis()

# 4) Produce all results and figures
# run_results()  # See results.R and plotting.R

# 5) Upload resuts to Acacia storage
run_acacia("download")  # See acacia.R

# Finish up
message("* Finished!")

