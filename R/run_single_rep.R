###########################################
#           SETONIX READY SCRIPT
###########################################


###########################################
#        SOURCE DEPENDENCIES
###########################################

source("dependencies.R")

###########################################
#        READ SLURM ARRAY ID
###########################################

args <- commandArgs(trailingOnly = TRUE)
task_id <- ifelse(length(args) > 0, as.numeric(args[1]), 1)

cat("Running SLURM task ID:", task_id, "\n")

###########################################
#        DEFINE PARAMETER GRID
###########################################

param_set <- expand.grid(
  lethal_effect    = c(TRUE, FALSE),
  complete_sterile = c(TRUE, FALSE)
) %>%
  mutate(scenario = row_number())

# Remove first combination if needed
param_set <- param_set[-1, ]

# Safety check
if (task_id > nrow(param_set)) {
  stop("SLURM task_id exceeds number of scenarios.")
}

params <- param_set[task_id, ]

###########################################
#        SET SCENARIO-SPECIFIC SEED
###########################################

set.seed(1000 + task_id)

###########################################
#        SETUP PARALLEL BACKEND
###########################################

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
cat("Using cores:", n_cores, "\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

###########################################
#        RUN REPLICATES IN PARALLEL
###########################################

results <- foreach(
  rep = 1:n_replicates,
  .packages = c("dplyr"),
  .combine = "c"
) %dopar% {
  
  # Independent replicate seed
  set.seed(100000 * task_id + rep)
  
  scenario_output <- run_model(
    n_patches         = patches,
    pop_patches       = pop_patches,
    n_per_patch       = n_per_patch,
    n_loci            = n_loci,
    init_frequency    = init_frequency,
    fecundity         = fecundity,
    carrying_capacity = carrying_capacity,
    decay             = decay,
    lambda            = lambda,
    lethal_effect     = params$lethal_effect,
    complete_sterile  = params$complete_sterile,
    linkage           = FALSE,
    sim_years         = sim_years,
    adjacency_matrix  = TRUE,
    dispersal_frac    = dispersal_frac
  )
  
  patch_stats <- scenario_output$patch_stats %>%
    mutate(
      scenario         = params$scenario,
      replicate        = rep,
      lethal_effect    = params$lethal_effect,
      complete_sterile = params$complete_sterile
    )
  
  genetic_stats <- scenario_output$genetic_data %>%
    mutate(
      scenario         = params$scenario,
      replicate        = rep,
      lethal_effect    = params$lethal_effect,
      complete_sterile = params$complete_sterile
    )
  
  list(
    patch   = patch_stats,
    genetic = genetic_stats
  )
}

stopCluster(cl)

cat("Parallel runs completed.\n")

###########################################
#        BIND OUTPUT
###########################################

all_patch_stats  <- bind_rows(lapply(results, `[[`, "patch"))
all_genetic_data <- bind_rows(lapply(results, `[[`, "genetic"))

###########################################
#        SAVE OUTPUT (SCENARIO-SPECIFIC)
###########################################

if (!dir.exists("output")) dir.create("output")

saveRDS(
  all_patch_stats,
  file = file.path("output",
                   paste0("patch_scenario_", task_id, ".rds"))
)

saveRDS(
  all_genetic_data,
  file = file.path("output",
                   paste0("genetic_scenario_", task_id, ".rds"))
)

cat("Output saved successfully for scenario", task_id, "\n")