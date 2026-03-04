

###########################################
#               PARAMETERS                #
###########################################

set.seed(230)

#This code sources all  required functions and sub models
source("dependencies.R")


###########################################
#            
#            RUN SINGLE SIMULATION      
#                   
###########################################
# Set working directory to sourced file




# -----------------------------
# Read environment variable from SLURM
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- ifelse(length(args) > 0, as.numeric(args[1]), 1)
cat("Running task ID:", task_id, "\n")



# -----------------------------------------------------------
# generating the parameter range manually with expand.grid
# ----------------------------------------------------------

param_set <- expand.grid(
  #dispersal_prob = c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05),
  lethal_effect = c(TRUE, FALSE),
  complete_sterile = c(TRUE, FALSE)
) %>%
  mutate(
    scenario = row_number()
  )


param_set <- param_set [-(1),]

# Get parameter set for this SLURM job
params <- param_set[task_id, ]


# -----------------------------------------
# Parallel backend setup (for replicates)
# -------------------------------------------

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)


# -----------------------------------------------
# Run model replicates in parallel
# -------------------------------------------

results <- foreach(rep = 1:n_replicates, .packages = c("dplyr")) %dopar% {
  scenario_output <- run_model(
    n_patches        = patches,
    pop_patches      = pop_patches,
    n_per_patch      = n_per_patch,
    n_loci           = n_loci,
    init_frequency   = init_frequency,
    fecundity        = fecundity,
    carrying_capacity= carrying_capacity,
    decay            = decay,
    lambda           = lambda,
    lethal_effect    = params$lethal_effect,
    complete_sterile = params$complete_sterile,
    linkage          = FALSE,
    sim_years        = sim_years,
    adjacency_matrix = TRUE,
    dispersal_frac   = dispersal_frac
  )
  
  
  # --- Add scenario + replicate details ---
  patch_stats <- scenario_output$patch_stats %>%
    mutate(
      scenario       = params$scenario,
      replicate      = rep,
      complete_sterile = params$complete_sterile,
      lethal_effect    = params$lethal_effect
    )
  
  genetic_stats <- scenario_output$genetic_data %>%
    mutate(
      scenario       = params$scenario,
      replicate      = rep,
      complete_sterile = params$complete_sterile,
      lethal_effect    = params$lethal_effect
    )
  
  #list(patch = patch_stats, genet = genetic_stats)
  # Append to collectors
  all_patch_stats <- append(all_patch_stats, list(patch_stats))
  all_genetic_data <- append(all_genetic_data, list(genetic_stats))
}

stopCluster(cl)


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

if (!dir.exists("output")) dir.create("output")
saveRDS(all_patch_stats, file = file.path("output", "step_data1000.rds"))
saveRDS(all_genetic_data, file = file.path("output", "step_genetic1000.rds"))

cat("Binding completed! Output saved", "\n")