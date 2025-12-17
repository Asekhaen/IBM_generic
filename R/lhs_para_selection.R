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


# --------------------------------------------------------------------------
# Using Latin Hypercube Sampling to generate parameter ranges to explore
# --------------------------------------------------------------------------

param_ranges <- data.frame(
  fecundity =  c(2, 10),
  n_loci          = c(100, 1000),
  init_frequency = c(0.01, 0.25),
  dispersal_prob  = c(0.001, 0.005)
)

set.seed(230)

lhs_sample <- randomLHS(n_samples, ncol(param_ranges))

param_set <- data.frame(
  fecundity = qinteger(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
  n_loci = qinteger(lhs_sample[,2], param_ranges[1,2], param_ranges[2,2]),
  init_frequency = qunif(lhs_sample[,3], param_ranges[1,3], param_ranges[2,3]),
  dispersal_prob = qunif(lhs_sample[,4], param_ranges[1,4], param_ranges[2,4]),
  scenario = 1:n_samples
)


# Get parameter set for this SLURM job
params <- param_set[task_id, ]

# -----------------------------
# Parallel back-end setup (for replicates)
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
    n_loci           = params$n_loci,
    init_frequency   = params$init_frequency,
    fecundity        = params$fecundity,
    carrying_capacity= carrying_capacity,
    lethal_effect    = FALSE,
    complete_sterile = FALSE,
    sim_years        = sim_years,
    dd_rate          = dd_rate,
    lambda           = lambda,
    adjacency_matrix = TRUE,
    dispersal_frac   = params$dispersal_frac,
    decay            = decay
  )
  
  patch_stats <- scenario_output$pop_stats |>
    mutate(
      scenario       = params$scenario,
      replicate      = rep,
      init_frequency = params$init_frequency,
      n_loci         = params$n_loci,
      fecundity      = params$fecundity,
      dispersal_frac = params$dispersal_frac
    )
  
  genetic_stats <- scenario_output$genetic_data |>
    mutate(
      scenario       = params$scenario,
      replicate      = rep,
      init_frequency = params$init_frequency,
      n_loci         = params$n_loci,
      fecundity      = params$fecundity,
      dispersal_frac = params$dispersal_frac
    )
  
  list(patch = patch_stats, genetics = genetic_stats)
}

stopCluster(cl)


# -----------------------------
# bind final and save outputs
# -----------------------------

all_patch_stats <- bind_rows(lapply(results, `[[`, "patch"))
all_genetic_stats <- bind_rows(lapply(results, `[[`, "genetic"))

if (!dir.exists("output")) dir.create("output")

saveRDS(all_patch_stats, file = sprintf("output/scenario_%03d_patch.rds", params$scenario))
saveRDS(all_genetic_stats, file = sprintf("output/scenario_%03d_genetics.rds", params$scenario))

cat("Scenario", params$scenario, "completed.\n")





# ================================================================
# combine scenario outputs from simulations into a single file
# ================================================================

# Make sure the "output" directory exists
if (!dir.exists("output")) stop("Output directory not found!")

# -----------------------------
#  Get file lists
# -----------------------------
patch_files <- list.files("output", pattern = "_patch\\.rds$", full.names = TRUE)
genetic_file <- list.files("output", pattern = "_genetic\\.rds$", full.names = TRUE)

cat("Found", length(patch_files), "patch files and", length(genetic_file), "genetic files.\n")

# -----------------------------
# Combine all patch results
# -----------------------------
all_patch_stats <- map_dfr(patch_files, readRDS)
cat("Combined patch data: ", nrow(all_patch_stats), " rows.\n")

# -----------------------------
# Combine all allele frequency results
# -----------------------------
all_genetic_stats <- map_dfr(genetic_file, readRDS)
cat("Combined allele data: ", nrow(all_genetic_stats), " rows.\n")

# -----------------------------
# Save final combined results
# -----------------------------
if (!dir.exists("combined")) dir.create("combined")

saveRDS(all_patch_stats, file = "combined/all_patch_stats_combined.rds")
saveRDS(all_genetic_stats, file = "combined/all_genetic_stats_combined.rds")

cat("✅ Combined results saved in 'combined/' directory.\n")

