
# load libraries needed


###########################################
#               PARAMETERS                #
###########################################

set.seed(230)


###########################################
#            
#            RUN SINGLE SIMULATION      
#                   
###########################################
# Set working directory to sourced file


 source("dependencies.R")



# -----------------------------
#  Single simulation 
# -----------------------------


# results <- run_model (
#   n_patches = patches,
#   pop_patches,
#   n_per_patch = n_per_patch,
#   n_loci = n_loci,
#   init_frequency = init_frequency,
#   fecundity = fecundity,
#   carrying_capacity = carrying_capacity,
#   prob_survival = prob_survival,
#   dd_rate = dd_rate,
#   decay = decay,
#   lambda = lambda,
#   lethal_effect = FALSE,
#   complete_sterile = TRUE,
#   sim_years = sim_years,
#   adjacency_matrix = TRUE,
#   dispersal_frac = dispersal_prob
# )




# ------------------------------------------------------
#  Running with varying parameters and replicates
# ---------------------------------------------------

# --------------------------------------------------------------------------
# Using Latin Hypercube Sampling to generate parameter ranges to explore
# --------------------------------------------------------------------------
# param_ranges <- data.frame(
#   # init_frequency = c(0.01, 0.30),
#   # n_loci          = c(50, 200),
#   dispersal_prob  = c(0.0001, 0.05)
# )
#
# lhs_sample <- randomLHS(n_samples, ncol(param_ranges))
# param_set <- data.frame(
#   # init_frequency = qunif(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
#   # n_loci = qinteger(lhs_sample[,2], param_ranges[1,2], param_ranges[2,2]),
#   dispersal_prob = qunif(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
#   scenario = 1:n_samples
# )

# -----------------------------------------------------------
# generating the parameter range manually with expand.grid
# ----------------------------------------------------------

param_set <- expand.grid(
  dispersal_prob = c(0.001, 0.005, 0.01, 0.025, 0.05),
  complete_sterile = c(TRUE, FALSE)
) %>%
  mutate(
    scenario = row_number()
  )

all_patch_stats <- list()
all_allele_frequency <- list()

for (i in 1:nrow(param_set)) {
  cat("Running parameter set:", param_set$scenario[i], "\n")
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
      complete_sterile = param_set$complete_sterile[i],
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
        dispersal_frac = param_set$dispersal_prob[i],
        complete_sterile = param_set$complete_sterile[i]
      )

    allele_stats <- scenario_output$allele_freq_per_locus |>
      mutate(
        #scenario       = param_set$scenario[i],
        replicate      = rep,
        #init_frequency = param_set$init_frequency[i],
        dispersal_frac = param_set$dispersal_prob[i],
        complete_sterile = param_set$complete_sterile[i]
      )

    # Append to collectors
    all_patch_stats <- append(all_patch_stats, list(patch_stats))
    all_allele_frequency <- append(all_allele_frequency, list(allele_stats))
  }
}

cat("Runs completed!", "\n")

# -----------------------------
# bind final outputs
# -----------------------------

cat("Binding output...", "\n")

all_patch_stats <- bind_rows(all_patch_stats)
all_allele_frequency <- bind_rows(all_allele_frequency)

# -----------------------------
# save bound outputs
# -----------------------------

if (!dir.exists("output")) dir.create("output")
saveRDS(all_patch_stats, file = file.path("output", "patch_data.rds"))
saveRDS(all_allele_frequency, file = file.path("output", "patch_freq.rds"))

cat("Binding completed! Output saved", "\n")
