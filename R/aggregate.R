############################################################
# AGGREGATE RESULTS
#
# Concatenate model output from all sims into single files.
#
############################################################

# ---------------------------------------------------------
# Aggregate results: history matching
# ---------------------------------------------------------
run_aggregate = function(sims) {
  
  message(" > Aggregating simulatons")

  # Remove redundant variables from sims info
  sims %<>% select(-run, -job_num)
  
  # Load model output for all simulations
  output = load_simulations(sims)

  # Concatenate seperate results for key model outcomes
  concatenate_result(sims, output, "patch_sizes")
  concatenate_result(sims, output, "allele_frequency")

  # TODO: Deal with 'final_pop'
}

# ---------------------------------------------------------
# Load all processed model output for a given set of sims
# ---------------------------------------------------------
load_simulations = function(sims) {
  
  # All postprocessed output files we wish to load
  output_files = paste0(o$pth$sims, sims$id, ".rds")
  
  # Load and concatenate into single datatable (in parallel)
  if (o$parallel_loading == TRUE)
    output = mclapply(output_files, load_attempt, mc.cores = o$cores)
  
  # Load and concatenate into single datatable (consecutively)
  if (o$parallel_loading == FALSE)
    output = lapply(output_files, load_attempt)
  
  # Name elements of list using sim IDs
  names(output) = sims$id
  
  return(output)
}

# ---------------------------------------------------------
# Attempt to load simulation result
# ---------------------------------------------------------
load_attempt = function(file) {
  
  # Initiate trivial output
  output = NULL
  
  # Load if file exists
  if (file.exists(file))
    output = readRDS(file)
  
  return(output)
}

# ---------------------------------------------------------
# Concatenate a given result into single datatable
# ---------------------------------------------------------
concatenate_result = function(sims, output, result) {

  # Extract outcomes for this result only
  result_dt = output %>%
      lapply(function(x) x[[result]]) %>%
      rbindlist(idcol = "id") %>%
      # Append simulation details...
      left_join(y  = sims, 
                by = "id") %>%
      relocate(names(sims), .before = 1)

  # Save aggregated, concatenated result
  save_rds(result_dt, "compiled", result)
}

