############################################################
# ANALYSIS
#
# Any stats analysis that operates on the concatenated
# simulations can be done here.
#
############################################################

# ---------------------------------------------------------
# Function for running analysis
# ---------------------------------------------------------
run_analysis = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Performing analysis")
  
  # ---- Load stuff ----
  
  # Simulation outputs: patch_sizes
  patch_sizes_dt = read_rds("compiled", "patch_sizes")
  
  # Simulation outputs: allele_frequency
  allele_frequency_dt = read_rds("compiled", "allele_frequency")
  
  # ---- Perform analysis here ----
  
  # TODO: Perform any stats analysis here...
  
}

