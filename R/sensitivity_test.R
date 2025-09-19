
# load libraries needed


###########################################
#               PARAMETERS                #
###########################################

set.seed(230)

# Source functions and parameters 
source("R/sub_functions.R")
source("R/generic_function.R")
source("R/parameters_generic.R")

packages <- c("ggplot2", 
              "dplyr",
              "tibble",
              "tidyr", 
              "readr", 
              "purrr",       # uses pmap to loop through different scenarios
              "furrr",       # multisession i.e. distribute work across many cores
              "progressr",
              "patchwork",
              "lhs",
              "gbm",
              "future.apply",
              "reshape2",
              "data.table")  


load_libraries(packages)


# -----------------------------
# 1. Parameter ranges for Latin Hypercube Sampling
# -----------------------------
param_ranges <- data.frame(
  init_frequency = c(0.1, 0.5),
  fecundity       = c(1, 10),
  n_loci          = c(5, 30),
  dispersal_prob  = c(0.0001, 0.01)
)

# -----------------------------
# 2. Generate and scale LHS Parameter Sets ---
# -----------------------------

n_samples <- 25
n_replicates <- 3
lhs_sample <- randomLHS(n_samples, ncol(param_ranges))
scale_lhs <- function(lhs, min, max) min + (max - min) * lhs
param_set <- data.frame(
  init_frequency = scale_lhs(lhs_sample[,1], param_ranges[1,1], param_ranges[2,1]),
  fecundity = round(scale_lhs(lhs_sample[,2], param_ranges[1,2], param_ranges[2,2])),
  n_loci = round(scale_lhs(lhs_sample[,3], param_ranges[1,3], param_ranges[2,3])),
  dispersal_prob = scale_lhs(lhs_sample[,4], param_ranges[1,4], param_ranges[2,4]),
  scenario = 1:n_samples
)




all_patch_stats <- list()
all_allele_frequency <- list()

for (i in 1:nrow(param_set)) {
  cat("Running sample", i, "/", n_samples, "\n")
  
  for (rep in 1:n_replicates) {
    
    scenario_output <- run_model(
      n_patches        = patches,
      pop_patches      = pop_patches, 
      n_per_patch      = n_per_patch,
      n_loci           = param_set$n_loci[i],
      init_frequency   = param_set$init_frequency[i],
      fecundity        = param_set$fecundity[i],
      carrying_capacity= carrying_capacity,
      lethal_effect    = FALSE,
      complete_sterile = FALSE,
      sim_years        = sim_years,
      prob_survival    = prob_survival,
      dd_rate          = dd_rate,
      overlapping      = TRUE,
      lambda           = lambda,
      adjacency_matrix = TRUE,
      dispersal_frac   = param_set$dispersal_prob[i],
      decay            = decay
    )
    
    # --- Add scenario + replicate details ---
    ps <- scenario_output$pop_stats |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        init_frequency = param_set$init_frequency[i],
        fecundity      = param_set$fecundity[i],
        n_loci         = param_set$n_loci[i],
        dispersal_prob = param_set$dispersal_prob[i]
      )
    
    af <- scenario_output$allele_freq_per_locus |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        init_frequency = param_set$init_frequency[i],
        fecundity      = param_set$fecundity[i],
        n_loci         = param_set$n_loci[i],
        dispersal_prob = param_set$dispersal_prob[i]
      )
    
    # Append to collectors
    all_patch_stats <- append(all_patch_stats, list(ps))
    all_allele_frequency <- append(all_allele_frequency, list(af))
  }
}

# -----------------------------
# Final bound outputs
# -----------------------------
all_patch_stats <- bind_rows(all_patch_stats)
all_allele_frequency <- bind_rows(all_allele_frequency)


#saveRDS(all_patch_stats, file = "C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\result_test.rds")


# aggregate replicates 
patch_stats <- all_patch_stats %>%
  group_by(patch, year, init_frequency, fecundity, n_loci, dispersal_prob) %>%
  summarise(
    pop_size = mean(pop_size, na.rm = TRUE),
    speed = mean(speed, na.rm = TRUE),
    freq = mean(freq, na.rm = TRUE),
    .groups = "drop"
  )


# fit a Boosted regression tree (BRT) using gbm package

# BRT for pop size
brt_pop <- gbm(
  formula = pop_size ~ init_frequency + fecundity + n_loci + dispersal_prob,
  data = patch_stats,
  distribution = "gaussian",  # use "bernoulli" for binary outcomes
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.5,
  train.fraction = 0.8,
  cv.folds = 5,
  verbose = FALSE
)

# shows relative influence of each predictor
summary(brt_pop) 

#BRT for allele frequency
brt_freq <- gbm(
  formula = freq ~ init_frequency + fecundity + n_loci + dispersal_prob,
  data = patch_stats,
  distribution = "gaussian",  # use "bernoulli" for binary outcomes
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.5,
  train.fraction = 0.8,
  cv.folds = 5,
  verbose = FALSE
)

# shows relative influence of each predictor
summary(brt_freq) 


# BRT Speed
brt_speed <- gbm(
  formula = speed ~ init_frequency + fecundity + n_loci + dispersal_prob,
  data = patch_stats,
  distribution = "gaussian",  # use "bernoulli" for binary outcomes
  n.trees = 5000,
  interaction.depth = 3,
  shrinkage = 0.01,
  bag.fraction = 0.5,
  train.fraction = 0.8,
  cv.folds = 5,
  verbose = FALSE
)

# shows relative influence of each predictor
summary(brt_speed) 


# plot(brt_pop, i.var = "init_frequency")
# plot(brt_pop, i.var = "fecundity")
# plot(brt_pop, i.var = "n_loci")
# plot(brt_pop, i.var = "dispersal_prob")



 





# results_list <- list()
# 
# for (i in 1:nrow(param_set)) {
#   cat("Running sample", i, "/", n_samples, "\n")
#   for (rep in 1:n_replicates) {
#     scenario_output <-run_model(
#       n_patches = patches,
#       pop_patches = pop_patches,
#       n_per_patch = n_per_patch,
#       n_loci = param_set$n_loci[i],
#       init_frequency = param_set$init_frequency[i],
#       fecundity = param_set$fecundity[i],
#       carrying_capacity = carrying_capacity,
#       lethal_effect = TRUE,
#       complete_sterile = FALSE,
#       sim_years = sim_years,
#       prob_survival = prob_survival,
#       dd_rate = dd_rate,
#       overlapping = TRUE,
#       lambda = lambda,
#       adjacency_matrix = TRUE,
#       dispersal_frac = param_set$dispersal_prob[i],
#       decay = decay
#     )
#     results_list[[length(results_list)+1]] <- c(param_set[i,], replicate=rep, scenario_output)
#   }
# }
# 
# pop_per_patch <- results_list[[1]]$patch_stats
# allele_per_patch <- results_list[[1]]$allele_freq_per_locus
# 


# all_patch_stats <- list()
# all_allele_frequency <- list()
# 
# for (i in 1:nrow(param_set)) {
#   cat("Running sample", i, "/", n_samples, "\n")
#   for (rep in 1:n_replicates) {
#     
#     scenario_output <- run_model(
#       n_patches = patches,
#       pop_patches = pop_patches, 
#       n_per_patch = n_per_patch,
#       n_loci = param_set$n_loci[i],
#       init_frequency = param_set$init_frequency[i],
#       fecundity = param_set$fecundity[i],
#       carrying_capacity = carrying_capacity,
#       lethal_effect = FALSE,
#       complete_sterile = FALSE,
#       sim_years = sim_years,
#       prob_survival = prob_survival,
#       dd_rate = dd_rate,
#       overlapping = TRUE,
#       lambda = lambda,
#       adjacency_matrix = TRUE,
#       dispersal_frac = param_set$dispersal_prob[i],
#       decay = decay
#     )
#     
#     # --- patch_stats table ---
#     all_patch_stats[[length(all_patch_stats) + 1]] <- scenario_output$patch_stats %>%
#       mutate(
#         scenario   = i,
#         replicate  = rep,
#         fecundity  = param_set$fecundity[i],
#         n_loci     = param_set$n_loci[i],
#         init_frequency = param_set$init_frequency[i],
#         dispersal_prob = param_set$dispersal_prob[i]
#       )
#     
#     # --- allele_frequency table ---
#     all_allele_frequency[[length(all_allele_frequency) + 1]] <- scenario_output$allele_frequency %>%
#       mutate(
#         scenario   = i,
#         replicate  = rep,
#         fecundity  = param_set$fecundity[i],
#         n_loci     = param_set$n_loci[i],
#         init_frequency = param_set$init_frequency[i],
#         dispersal_prob = param_set$dispersal_prob[i]
#       )
#   }
# }
# 
# # Combine all into single dataframes
#  
# all_patch_stats <- bind_rows(all_patch_stats)
# all_allele_frequency <- bind_rows(all_allele_frequency)
 


