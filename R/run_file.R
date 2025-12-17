
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
#  Single run 
# -----------------------------


results <- run_model (
  patches = patches,
  pop_patches,
  n_per_patch = n_per_patch,
  n_loci = n_loci,
  init_frequency = init_frequency,
  fecundity = fecundity,
  carrying_capacity = carrying_capacity,
  decay = decay,
  lambda = lambda,
  lethal_effect = FALSE,
  complete_sterile = FALSE,
  linkage = FALSE,
  sim_years = sim_years,
  adjacency_matrix = FALSE,
  dispersal_frac = dispersal_frac
)




# ------------------------------------------------------
#  multiple runs, varying parameters and replicates
# ---------------------------------------------------

# -----------------------------------------------------------
# generating the parameter range manually with expand.grid
# ----------------------------------------------------------

# source("dependencies.R")
# 
# 
# param_set <- expand.grid(
#   #dispersal_prob = c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05),
#   complete_sterile = c(TRUE, FALSE)
# ) %>%
#   mutate(
#     scenario = row_number()
#   )
# 
# all_patch_stats <- list()
# all_genetic_data <- list()
# 
# for (i in 1:nrow(param_set)) {
#   cat("Running parameter set:", param_set$scenario[i], "\n")
#   for (rep in 1:n_replicates) {
# 
#     scenario_output <- run_model(
#       patches        = patches,
#       #pop_patches      = pop_patches,
#       n_per_patch      = n_per_patch,
#       n_loci           = n_loci,
#       init_frequency   = init_frequency,
#       fecundity        = fecundity,
#       carrying_capacity= carrying_capacity,
#       lethal_effect    = FALSE,
#       complete_sterile = param_set$complete_sterile[i],
#       linkage = FALSE,
#       sim_years        = sim_years,
#       lambda           = lambda,
#       adjacency_matrix = TRUE,
#       dispersal_frac   = dispersal_frac,
#       decay            = decay
#     )
# 
#     # --- Add scenario + replicate details ---
#     patch_stats <- scenario_output$patch_stats |>
#       mutate(
#         #scenario       = param_set$scenario[i],
#         replicate      = rep,
#         #init_frequency = param_set$init_frequency[i],
#         #dispersal_frac = dispersal_frac,
#         complete_sterile = param_set$complete_sterile[i]
#       )
# 
#     genetic_stats <- scenario_output$genetic_data |>
#       mutate(
#         #scenario       = param_set$scenario[i],
#         replicate      = rep,
#         #init_frequency = param_set$init_frequency[i],
#         #dispersal_frac = dispersal_frac,
#         complete_sterile = param_set$complete_sterile[i]
#       )
# 
#     # Append to collectors
#     all_patch_stats <- append(all_patch_stats, list(patch_stats))
#     all_genetic_data <- append(all_genetic_data, list(genetic_stats))
#   }
# }
# 
# cat("Runs completed!", "\n")
# 
# # -----------------------------
# # bind final outputs
# # -----------------------------
# 
# cat("Binding output...", "\n")
# 
# all_patch_stats <- bind_rows(all_patch_stats)
# all_genetic_data <- bind_rows(all_genetic_data)
# 
# # -----------------------------
# # save bound outputs
# # -----------------------------
# 
# if (!dir.exists("output")) dir.create("output")
# saveRDS(all_patch_stats, file = file.path("output", "patch_data.rds"))
# saveRDS(all_genetic_data, file = file.path("output", "egentic_data.rds"))
# 
# cat("Binding completed! Output saved", "\n")
