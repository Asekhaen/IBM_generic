
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


# Get parameter set for this SLURM job
params <- param_set[task_id, ]


# -----------------------------------------
# Parallel backend setup (for replicates)
# -------------------------------------------

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cl <- makeCluster(n_cores, type = "PSOCK")
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
      linkage = TRUE,
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
    
    list(patch = patch_stats, genet = genetic_stats)
  }

stopCluster(cl)

cat("Runs completed!", "\n")


# -----------------------------
# bind and save final outputs
# -----------------------------

all_patch_stats <- bind_rows(lapply(results, `[[`, "patch"))
all_genetic_stats <- bind_rows(lapply(results, `[[`, "genet"))

if (!dir.exists("output")) dir.create("output")

saveRDS(all_patch_stats, file = sprintf("output/scenario_%03d_patch.rds", params$scenario))
saveRDS(all_genetic_stats, file = sprintf("output/scenario_%03d_allele.rds", params$scenario))
