###########################################################
# OPTIONS
#
# Define a whole load of simulation settings and assumptions.
#
############################################################

# ---------------------------------------------------------
# Set a series of options for the analysis
# ---------------------------------------------------------
set_options = function(do_step = NA, quiet = FALSE) {
  
  if (quiet == FALSE) message("* Setting options")
  
  # ---- Basics ----
  
  # Reset R's most annoying default options
  # default_R_options()
  
  # Details of the user and their setup
  o = identify_user()
  
  # Append which step(s) we're running
  o$do_step = do_step
  
  # Output version (linked to output directory name)
  o$output_version = "v01"
  
  # Set directories
  o = set_dirs(o)  # See directories.R
  
  # ---- Analysis settings ----
  
  # Things that vary in simulations...
  
  # Dispersal probability
  o$dispersal_prob <- c(
    high = 0.01,
    med  = 0.001,
    low  = 0.0005)
  
  # Initial frequency
  o$init_frequency <- c(
    high = 0.25, 
    low  = 0.1)
  
  # xxx....
  o$lethal_effect = c(
    yes = TRUE, 
    no  = FALSE)
  
  # Number of offspring per day per female mosquito
  o$fecundity <- c(
    default = 5)
  
  # ---- Simulation settings ----
  
  # Things that remain fixed in each simulation
  
  # Model options
  
  o$prob_survival <- 0.7
  o$dd_rate <- 0.0001
  o$patches <- 7                               # Number of patches
  o$n_per_patch <- c(10000,0,0,0,0,0,0)    # Initial number of individuals per patch
  o$beta <- 100                           # the adult male population size at which the daily probability of mating is 0.5.
  o$sim_years <-20                         # Number of simulation in days
  o$carrying_capacity = 10000             # carrying capacity  
  
  # dispersal parameters
  o$lambda <- 0.1
    
  # Genetic (load) & drive parameters
  o$n_loci <- 10
  o$decay <- 0.5  
  
  # Set flags for certain model features
  o$complete_sterile = TRUE
  o$overlapping = TRUE
  
  # ---- Seasonality settings ----
  
  # Seasonal settings
  o$season = c("seasonal", "perennial")
  
  # ---- Time settings ----
  
  # Number of years over which to assess results
  o$years_assess = 15
  
  # Years of analysis
  o$year_start   = 2000
  o$year_current = 2024
  o$year_final   = 2040
  
  # ---- Uncertainty settings ----
  
  # Number of seeds per scenario
  o$n_seeds = 5
  
  # Results summary prediction interval
  o$pred_int = 95  # As percentage
  
  # ---- Post-processing settings ----
  
  # Load all results in parallel
  o$parallel_loading = TRUE
  
  # ---- Cluster settings ----
  
  # Flag for overwriting any existing simulations
  o$overwrite = FALSE
  
  # Number of times to re-run failed simulations before aggregating
  o$rerun_sims = 3
  
  # HPC-specific settings and options
  o$job_settings = list(
    # For running on Setonix...
    setonix = list(
      account   = o$hpc$project,
      partition = "work",
      export    = "NONE"),
    # For running on SciCore...
    scicore = list(
      account   = "penny",
      qos       = "30min"))
  
  # Job size settings (HPC agnostic)
  o$job_size = list(
    'nodes'         = 1,
    'ntasks'        = 128, 
    'cpus-per-task' = 1,
    'mem-per-cpu'   = "500MB", 
    'time'          = "00:20:00")
  
  # Set an upper limit for array jobs that can be run at any one time
  o$job_limit = 1000
  
  # Define names for cluster log and error files
  o$log_file = "cluster_log.txt"
  o$err_file = "cluster_error.txt"
  
  # Action to take if user is already running cluster jobs
  o$cluster_conflict_action = "error"  # Set to 'none' to turn off
  
  # ---- Results and figure flags ----
  
  # Plotting flags (toggle on or off)
  o$plot = list(
    setting     = T,  # Setting-specific details
    samples     = T,  # Parameter sampling visuals
    eir2prev    = T,  # EIR-prevalence fit results
    projections = F,  # Burden over time
    impact      = T)  # Public health impact results
  
  # ---- Plotting options ----
  
  # Saved figure size
  o$save_width  = 9
  o$save_height = 6
  
  # Units of figures sizes
  o$save_units = "in"
  
  # Plotting resolution (in dpi)
  o$save_resolution = 300
  
  # Image format for saving figure
  # 
  # NOTE: Use a character vector to save with multiple formats at once
  o$figure_format = "png" # Classic options: "png", "pdf", or "svg"
  
  # ---- Retain options ----
  
  # Save options to output dir for future reference
  saveRDS(o, paste0(o$pth$output, "options.rds"))
  
  # Report back some key options
  message(" > Analysis version: ", o$output_version)
  
  return(o)
}

