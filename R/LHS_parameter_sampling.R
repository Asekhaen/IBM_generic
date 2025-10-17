
# load libraries needed


###########################################
#               PARAMETERS                #
###########################################

#set.seed(230)

# Source functions, parameters, etc. from dependencies.R file
source("dependencies.R")

# 
# packages <- c("ggplot2", 
#               "dplyr",
#               "tibble",
#               "tidyr", 
#               "readr", 
#               "purrr",       # uses pmap to loop through different scenarios
#               "furrr",       # multisession i.e. distribute work across many cores
#               "progressr",
#               "patchwork",
#               "lhs",
#               "gbm",
#               "future.apply",
#               "reshape2",
#               "data.table")  
# 
# 
# load_libraries(packages)



# -----------------------------
# 1. Parameter ranges for Latin Hypercube Sampling
# -----------------------------
param_ranges <- data.frame(
  init_frequency = c(0.01, 0.25),
  n_loci          = c(100, 1000),
  fecundity       = c(2, 10),
  dispersal_prob  = c(0.001, 0.05)
)


# -----------------------------
# 2. Generate and scale LHS Parameter Sets ---
# -----------------------------
# sensitivity analysis


lhs_sample <- randomLHS(n_samples, ncol(param_ranges))
param_set <- data.frame(
  init_frequency = qunif(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
  n_loci = qinteger(lhs_sample[,2], param_ranges[1,2], param_ranges[2,2]),
  fecundity = qinteger(lhs_sample[,3], param_ranges[1,3], param_ranges[2,3]),
  dispersal_prob = qunif(lhs_sample[,4], param_ranges[1,4], param_ranges[2,4]),
  scenario = 1:n_samples
)


all_patch_stats <- list()
all_allele_frequency <- list()

for (i in 1:nrow(param_set)) {
  cat("Running sample", i, "/", n_samples, "\n")
  for (rep in 1:n_replicates) {
    
    scenario_output <- run_model(
      n_patches        = patches,
      #pop_patches      = pop_patches, 
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
        scenario       = param_set$scenario[i],
        replicate      = rep,
        init_frequency = param_set$init_frequency[i],
        n_loci         = param_set$n_loci[i],
        fecundity      = param_set$fecundity[i],
        dispersal_prob = param_set$dispersal_prob[i]
      )
    
    allele_stats <- scenario_output$allele_freq_per_locus |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        init_frequency = param_set$init_frequency[i],
        n_loci         = param_set$n_loci[i],
        fecundity      = param_set$fecundity[i],
        dispersal_prob = param_set$dispersal_prob[i]
      )
    
    # Append to collectors
    all_patch_stats <- append(all_patch_stats, list(patch_stats))
    all_allele_frequency <- append(all_allele_frequency, list(allele_stats))
  }
}


# -----------------------------
# bind final outputs
# -----------------------------

all_patch_stats <- bind_rows(all_patch_stats)
all_allele_frequency <- bind_rows(all_allele_frequency)


if (!dir.exists("output")) dir.create("output")

saveRDS(all_patch_stats, file = file.path("output", "parameter.rds"))
saveRDS(all_allele_frequency, file = file.path("output", "parameter_freq.rds"))

