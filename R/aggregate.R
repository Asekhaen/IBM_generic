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
  
  # Load model output for all simulations
  output_list = load_simulations(sims)
  
  browser()
  
  # Extract and concatenate patch size results
  patch_sizes_dt = output_list %>%
    lapply(function(x) x$patch_sizes) %>%
    rbindlist(idcol = "id")

  # Extract and concatenate allele frequency results
  allele_frequency_dt = output_list %>%
    lapply(function(x) x$allele_frequency) %>%
    rbindlist(idcol = "id")
  
  # Save aggregated, concatenated results
  save_rds(patch_sizes_dt,      "compiled", "patch_sizes")
  save_rds(allele_frequency_dt, "compiled", "allele_frequency")
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

