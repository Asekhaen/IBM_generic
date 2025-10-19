# ===============================
# run_parallel.R
# ===============================

# source dependencies

source("dependencies.R")

# -----------------------------
# 1. Read environment variable from SLURM
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- ifelse(length(args) > 0, as.numeric(args[1]), 1)
cat("Running task ID:", task_id, "\n")

# -----------------------------
# 2. Parameter ranges
# -----------------------------
param_ranges <- data.frame(
  init_frequency = c(0.01, 0.25),
  n_loci         = c(100, 1000),
  fecundity      = c(2, 10),
  dispersal_prob = c(0.001, 0.05)
)


# -----------------------------
# 3. Latin Hypercube Sampling (predefined parameter sets)
# -----------------------------
set.seed(123)

lhs_sample <- randomLHS(n_samples, ncol(param_ranges))
param_set <- data.frame(
  init_frequency = qunif(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
  n_loci         = qinteger(lhs_sample[,2], param_ranges[1,2], param_ranges[2,2]),
  fecundity      = qinteger(lhs_sample[,3], param_ranges[1,3], param_ranges[2,3]),
  dispersal_prob = qunif(lhs_sample[,4], param_ranges[1,4], param_ranges[2,4]),
  scenario       = 1:n_samples
)

# Get parameter set for this SLURM job
params <- param_set[task_id, ]

# -----------------------------
# 4. Parallel backend setup (for replicates)
# -----------------------------

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)



# -----------------------------
# 5. Run model replicates in parallel
# -----------------------------
results <- foreach(rep = 1:n_replicates, .packages = c("dplyr")) %dopar% {
  
  scenario_output <- run_model(
    n_patches        = patches,
    n_per_patch      = n_per_patch,
    n_loci           = params$n_loci,
    init_frequency   = params$init_frequency,
    fecundity        = params$fecundity,
    carrying_capacity= carrying_capacity,
    lethal_effect    = FALSE,
    complete_sterile = TRUE,
    sim_years        = sim_years,
    dd_rate          = dd_rate,
    lambda           = lambda,
    adjacency_matrix = TRUE,
    dispersal_frac   = params$dispersal_prob,
    decay            = decay
  )
  
  patch_stats <- scenario_output$pop_stats |>
    mutate(
      scenario       = params$scenario,
      replicate      = rep,
      init_frequency = params$init_frequency,
      n_loci         = params$n_loci,
      fecundity      = params$fecundity,
      dispersal_prob = params$dispersal_prob
    )
  
  allele_stats <- scenario_output$allele_freq_per_locus |>
    mutate(
      scenario       = params$scenario,
      replicate      = rep,
      init_frequency = params$init_frequency,
      n_loci         = params$n_loci,
      fecundity      = params$fecundity,
      dispersal_prob = params$dispersal_prob
    )
  
  list(patch = patch_stats, allele = allele_stats)
}

stopCluster(cl)

# -----------------------------
# 6. Bind and save results
# -----------------------------
all_patch_stats <- bind_rows(lapply(results, `[[`, "patch"))
all_allele_frequency <- bind_rows(lapply(results, `[[`, "allele"))

if (!dir.exists("output")) dir.create("output")

saveRDS(all_patch_stats, file = sprintf("output/scenario_%03d_patch.rds", params$scenario))
saveRDS(all_allele_frequency, file = sprintf("output/scenario_%03d_allele.rds", params$scenario))

cat("Scenario", params$scenario, "completed.\n")
