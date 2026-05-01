
# load libraries needed


#set.seed(230)


###########################################
#      RUNNING SINGLE SIMULATION
###########################################

# load source file. This include the main and sub models, and the parameter values

source("R/dependencies.R")



# -----------------------------
#  Single simulation. 
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
  complete_sterile = TRUE,
  linkage = FALSE,
  sim_years = sim_years,
  adjacency_matrix = TRUE,
  dispersal_frac = dispersal_frac
)




#The results is a list that tracks the population dynamics and the genetic stats



# ------------------------------------------------------
#  multiple runs, varying parameters and replicates
# ---------------------------------------------------


source("R/dependencies.R")


param_set <- expand.grid(
  dispersal_frac = c(0.001, 0.0025, 0.005, 0.01),
  # n_load = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5),
  # n_loci = c(1, 10, 100, 1000),
  lethal_effect = c(TRUE, FALSE),
  complete_sterile = c(TRUE, FALSE)
) |>
  mutate(
    scenario = row_number(),
    init_freq = calc_q(n_load,n_loci)
  )


#this line of code is used to remove lethal_effect = TRUE and complete_sterile = TRUE
#can be modified

param_set <- param_set [-(1:4),]


#create folder to save outputs if it doesn't exist

if (!dir.exists("R/output")) dir.create("R/output")
write_csv(param_set, file = "R/output/param_two_patch.csv")

#run multiple simulations
all_patch_stats <- list()
all_genetic_data <- list()

for (i in 1:nrow(param_set)) {

  cat("Running parameter set:", param_set$scenario[i], "\n")
  
  for (rep in 1:n_replicates) {

    scenario_output <- run_model (
      patches = patches,
      pop_patches,
      n_per_patch = n_per_patch,
      n_loci = n_loci, # param_set$n_loci[i],
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
      dispersal_frac = param_set$dispersal_frac[i]
    )

    # --- Add scenario + replicate details ---
    patch_stats <- scenario_output$patch_stats |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        lethal_effect = param_set$lethal_effect[i],
        complete_sterile = param_set$complete_sterile[i],
        #n_loci = param_set$n_loci[i],
        #n_load = param_set$n_load[i],
        #init_frequency = param_set$init_freq[i]
        dispersal_frac = param_set$dispersal_frac[i]
      )

    genetic_stats <- scenario_output$genetic_data |>
      mutate(
        scenario       = param_set$scenario[i],
        replicate      = rep,
        lethal_effect = param_set$lethal_effect[i],
        complete_sterile = param_set$complete_sterile[i],
        #n_loci = param_set$n_loci[i],
        # n_load = param_set$n_load[i],
        # init_frequency = param_set$init_freq[i],
        dispersal_frac = param_set$dispersal_frac[i]
      )

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
# save bound outputs. This saves the population dynamics = all_patch_stats 
# and genetics data = all_genetic_data
# -----------------------------

if (!dir.exists("R/output")) dir.create("R/output")
saveRDS(all_patch_stats, file = file.path("R/output", "step_data.rds"))
saveRDS(all_genetic_data, file = file.path("R/output", "step_genetic.rds"))

cat("Binding completed! Output saved", "\n")



