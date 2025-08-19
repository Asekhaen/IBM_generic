############################################################
# TEST
#
# Run a single simulation to test the model.
#
############################################################

# ---------------------------------------------------------
# Parent function for creating and running all simulations
# ---------------------------------------------------------
run_test = function(idx = 1) {
  
  # Only continue if specified by do_step
  if (!is.element(0, o$do_step)) return()
  
  message("* Running single test simulation")
  
  # Details of all simulations
  all_sims = get_simulations(test = TRUE)
  
  # Select simulation assocaited with idx
  sim = as.list(all_sims[idx, ])
  
  # Run model for this simulation
  result = run_model(sim)  # See generic_function.R
  
  message(" > Simulation complete")
}

