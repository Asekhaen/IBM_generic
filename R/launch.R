############################################################
# LAUNCH
#
# Launch pipeline for mosquito individual-based model.
#
############################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

message("Running mosquito IBM")

# Load options into global environment
o = set_options(do_step = 4)  # See options.R

# 1) Run simulations
run_simulations()  # See simulations.R

# 2) Perform stats analysis
# run_analysis()

# 3) Produce all results and figures
# run_results()  # See results.R and plotting.R

# 4) Upload resuts to Acacia storage
run_acacia("upload")  # See acacia.R

# Finish up
message("* Finished!")

