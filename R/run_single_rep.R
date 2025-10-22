# ===============================
# run_parallel.R
# ===============================

# source dependencies

source("dependencies.R")


# -----------------------------
# 4. Single with replicate setup (for replicates)
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
    n_loci           = pn_loci,
    init_frequency   = init_frequency,
    fecundity        = fecundity,
    carrying_capacity= carrying_capacity,
    lethal_effect    = FALSE,
    complete_sterile = TRUE,
    sim_years        = sim_years,
    dd_rate          = dd_rate,
    lambda           = lambda,
    adjacency_matrix = TRUE,
    dispersal_frac   = dispersal_prob,
    decay            = decay
  )
  
  patch_stats <- scenario_output$pop_stats |> mutate(replicate = rep)
  allele_stats <- scenario_output$allele_freq_per_locus |> mutate(replicate = rep)

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

