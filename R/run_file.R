###########################################
#            
#      Running different dispersal rate 
#                   
###########################################

source("dependencies.R")

# Detect number of CPUs from SLURM
args <- commandArgs(trailingOnly = TRUE)
task_id <- ifelse(length(args) > 0, as.numeric(args[1]), 1)
cat("Running task ID:", task_id, "\n")



# -----------------------------
# Parameter ranges for Latin Hypercube Sampling
# -----------------------------
param_ranges <- data.frame(
  # init_frequency = c(0.01, 0.30),
  # n_loci          = c(50, 200),
  dispersal_prob  = c(0.0009, 0.05)
)



# -----------------------------
# Generate and scale LHS Parameter Sets ---
# -----------------------------


lhs_sample <- randomLHS(n_samples, ncol(param_ranges))
param_set <- data.frame(
  # init_frequency = qunif(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
  # n_loci = qinteger(lhs_sample[,2], param_ranges[1,2], param_ranges[2,2]),
  dispersal_prob = qunif(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
  scenario = 1:n_samples
)


# Get parameter set for this SLURM job
params <- param_set[task_id, ]


# -----------------------------
# Parallel backend setup (for replicates)
# -----------------------------

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)


# -----------------------------
# Run model replicates in parallel
# ----------------------------- 

results <- foreach(rep = 1:n_replicates, .packages = c("dplyr")) %dopar% {
    scenario_output <- run_model(
      n_patches        = patches,
      n_per_patch      = n_per_patch,
      n_loci           = n_loci,
      init_frequency   = init_frequency,
      fecundity        = fecundity,
      carrying_capacity= carrying_capacity,
      lethal_effect    = FALSE,
      complete_sterile = FALSE,
      sim_years        = sim_years,
      dd_rate          = dd_rate,
      lambda           = lambda,
      adjacency_matrix = TRUE,
      dispersal_frac   = param_set$dispersal_prob[i],
      decay            = decay
    )
    
    # --- Add scenario + replicate details ---
    patch_stats <- scenario_output$pop_stats |>
      mutate(
        #scenario       = param_set$scenario[i],
        replicate      = rep,
        #init_frequency = param_set$init_frequency[i],
        #n_loci         = param_set$n_loci[i],
        dispersal_prob = param_set$dispersal_prob[i]
      )
    
    allele_stats <- scenario_output$allele_freq_per_locus |>
      mutate(
        #scenario       = param_set$scenario[i],
        replicate      = rep,
        #init_frequency = param_set$init_frequency[i],
        #n_loci         = param_set$n_loci[i],
        dispersal_prob = param_set$dispersal_prob[i]
      )
    
    list(patch = patch_stats, allele = allele_stats)
}


stopCluster(cl)

# -----------------------------
# Bind and save results
# -----------------------------
all_patch_stats <- bind_rows(lapply(results, `[[`, "patch"))
all_allele_frequency <- bind_rows(lapply(results, `[[`, "allele"))

if (!dir.exists("output")) dir.create("output")

saveRDS(all_patch_stats, file = sprintf("output/scenario_%03d_patch.rds", params$scenario))
saveRDS(all_allele_frequency, file = sprintf("output/scenario_%03d_allele.rds", params$scenario))

cat("Scenario", params$scenario, "completed.\n")
