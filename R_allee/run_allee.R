
# load libraries needed


###########################################
#               PARAMETERS                #
###########################################

#set.seed(230)


###########################################
#            
#            RUN SINGLE SIMULATION      
#                   
###########################################
# Set working directory to sourced file


 
# source("R_allee/dependencies_allee.R")
# 
# 
# 
# # -----------------------------
# #  Single run
# # -----------------------------
# 
# 
# results <- run_model (
#   patches = patches,
#   pop_patches,
#   n_per_patch = n_per_patch,
#   n_loci = n_loci,
#   init_frequency = init_frequency,
#   fecundity = fecundity,
#   carrying_capacity = carrying_capacity,
#   lambda = lambda,
#   lethal_effect = FALSE,
#   complete_sterile = TRUE,
#   linkage = FALSE,
#   sim_years = sim_years
# )



# ------------------------------------------------------
#  multiple runs, varying parameters and replicates
# ---------------------------------------------------

# -----------------------------------------------------------
# generating the parameter range manually with expand.grid
# ----------------------------------------------------------

source("R_allee/dependencies_allee.R")


param_set <- expand.grid(
  n_load = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5),
  n_loci = c(1, 10, 100, 1000),
  lethal_effect = c(TRUE, FALSE),
  complete_sterile = c(TRUE, FALSE)
) |>
  mutate(
    scenario = row_number(),
    init_freq = calc_q(n_load,n_loci)
  )

param_set <- param_set [-(1:24),]

write_csv(param_set, file = "R_allee/output/param_set_allee.csv")

all_patch_stats <- list()
all_genetic_data <- list()

for (i in 1:nrow(param_set)) {
  cat("Running parameter set:", param_set$scenario[i], "\n")
  
  # cal init freq for each scenario
  # init_frequency <- calc_q(param_set$n_load[i], param_set$n_loci[i])
  
  for (rep in 1:n_replicates) {
    
    scenario_output <- run_model (
      patches = patches,
      # pop_patches,
      n_per_patch = n_per_patch,
      n_loci = param_set$n_loci[i],
      init_frequency = param_set$init_freq[i],
      fecundity = fecundity,
      carrying_capacity = carrying_capacity,
      #decay = decay,
      lambda = lambda,
      lethal_effect = param_set$lethal_effect[i],
      complete_sterile = param_set$complete_sterile[i],
      linkage = FALSE,
      sim_years = sim_years,
      adjacency_matrix = TRUE,
    )
    
    # --- Add scenario + replicate details ---
    patch_stats <- scenario_output$patch_stats |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        lethal_effect = param_set$lethal_effect[i],
        complete_sterile = param_set$complete_sterile[i],
        n_loci = param_set$n_loci[i],
        n_load = param_set$n_load[i],
        init_frequency = param_set$init_freq[i],
        
      )
    
    genetic_stats <- scenario_output$genetic_data |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        lethal_effect = param_set$lethal_effect[i],
        complete_sterile = param_set$complete_sterile[i],
        n_loci = param_set$n_loci[i],
        n_load = param_set$n_load[i],
        init_frequency = param_set$init_freq[i],
      )
    
    # Append to collectors
    all_patch_stats <- append(all_patch_stats, list(patch_stats))
    all_genetic_data <- append(all_genetic_data, list(genetic_stats))
  }
}

cat("Runs completed!", "\n")

# -----------------------------
# bind final outputs
# -----------------------------

cat("Binding output...", "\n")

all_patch_stats <- bind_rows(all_patch_stats)
all_genetic_data <- bind_rows(all_genetic_data)

# -----------------------------
# save bound outputs
# -----------------------------

if (!dir.exists("R_allee/output")) dir.create("R_allee/output")
saveRDS(all_patch_stats, file = file.path("R_allee/output", "allee_effect.rds"))
saveRDS(all_genetic_data, file = file.path("R_allee/output", "allee_genetic.rds"))

cat("Binding completed! Output saved", "\n")



# two_patch_data <- split(two_patch_data, two_patch_data$scenario)
# two_genetic_data <- split(two_genetic, two_genetic$scenario)
#
# if (!dir.exists("step")) dir.create("step")
# saveRDS(step_stone_param, file ="output/step/step_stone_param.rds" )
# saveRDS(two_patch_data, file ="step/two_patch_data.rds" )
# saveRDS(two_genetic_data, file ="step/two_genetic_data.rds" )

