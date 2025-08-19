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
o = set_options(do_step = 1 : 4)  # See options.R

# 0) Run single test simualtion
run_test()  # See test.R

# 1) Run simulations
run_simulations()  # See simulations.R

# 2) Perform stats analysis
run_analysis()  # See analysis.R

# 3) Produce figures
run_plotting()  # See plotting.R

# 4) Upload resuts to Acacia storage
run_acacia("upload")  # See acacia.R

# Finish up
message("* Finished!")

