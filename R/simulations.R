############################################################
# SIMULATIONS
#
# Create, simulate, and post-process all OM simulations.
#
############################################################

# ---------------------------------------------------------
# Parent function for creating and running all simulations
# ---------------------------------------------------------
run_simulations = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Running simulations")
  
  # Repeat this step several times
  for (i in 0 : o$rerun_sims + 1) {
    
    message("* Cluster submission attempt: ", i)
    
    # Generate full list of simulations to run
    sims = get_simulations()
    
    # Number of jobs to be run
    n_jobs = sum(sims$run)
    
    # We can break out if no jobs to run
    if (n_jobs == 0)
      break()
    
    # Submit all jobs to the cluster (see auxiliary.R)
    submit_cluster_jobs(n_jobs, "submit.sh", "job_simulation")
    
    # Throw an error if any cluster jobs failed (see auxiliary.R)
    stop_if_errors(o$pth$log, o$err_file, err_tol = 0)
  }
  
  # Concatenate outputparam_sets
  run_aggregate(sims) # See aggregate.R
}

# ---------------------------------------------------------
# Generate full list of simulations to run based on analysis phase
# ---------------------------------------------------------
get_simulations = function() {
  
  # Construct simulation set
  sims = expand_grid(
    dispersal_prob = names(o$dispersal_prob), 
    init_frequency = names(o$init_frequency), 
    lethal_effect  = names(o$lethal_effect),
    fecundity      = names(o$fecundity),
    seed = 1 : o$n_seeds)
  
  # Create IDs for each sim
  ids = unite(sims, "ids", everything(), sep = "_")$ids
  
  # Append machine-readable simulation IDs
  sims %<>% 
    mutate(id = ids, .before = 1)
  
  # ---- Skip existing sims ----
  
  message(" > Identifying previously completed simulations")
  
  # Extract IDs of sims that have already been run
  exist_id = str_remove(list.files(o$pth$sims), ".rds$")
  
  # Logical whether simulation should be run / rerun
  run_sim = !(sims$id %in% exist_id)
  if (o$overwrite) run_sim[] = TRUE
  
  # Skip any existing sims (unless overwriting)
  sims %<>%
    cbind(run = run_sim) %>%
    mutate(job_num = cumsum(run * 1), 
           job_num = ifelse(run, job_num, NA))
  
  # Save scenario dataframe to file
  save_rds(sims, "compiled", "simulations")
  
  # ---- Display number of sims ----
  
  # Number of sims
  n_total = nrow(sims)
  n_run   = sum(sims$run)
  
  # Report total number of sims
  message(" > Total number of simulations: ", thou_sep(n_total))
  
  # Report number of sims we'll run now
  message("  - Skipping: ", thou_sep(n_total - n_run))
  message("  - Simulating: ", thou_sep(n_run))
  
  return(sims)
}

# ---------------------------------------------------------
# All steps to actually simulate the model
# ---------------------------------------------------------
job_simulation = function(job_id) {
  
  # Load details of all simulations
  all_sims = read_rds("compiled", "simulations")
  
  # Select simulation assocaited with job_id
  sim = all_sims %>% 
    filter(job_num == !!job_id) %>%
    select(-job_num, -run) %>%
    as.list()
  
  message(" > Simulating: ", sim$id)
  
  # Run model for this simulation
  run_model(sim)  # See generic_function.R
  
  message(" > Simulation complete")
}

