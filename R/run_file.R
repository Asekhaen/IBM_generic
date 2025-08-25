
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
              "data.table")   # shows the progress


load_libraries(packages)


###########################################
#            
#            RUN SINGLE SIMULATION      
#                   
###########################################

# results <- run_model (n_patches = patches,
#            pop_patches,
#            n_per_patch = n_per_patch,
#            n_loci = n_loci,
#            init_frequency = init_frequency,
#            fecundity = fecundity,
#            carrying_capacity = carrying_capacity,
#            prob_survival = prob_survival,
#            decay = decay,
#            lambda = lambda,
#            lethal_effect = FALSE,
#            complete_sterile = TRUE,
#            sim_years = sim_years,
#            overlapping = TRUE,
#            adjacency_matrix = TRUE,
#            dispersal_frac = dispersal_prob)
# 
# if (!dir.exists("output")) dir.create("output")
# saveRDS(output, file = "C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\results.rds")


  
  
  ###################################################
  #
  #   RUN PARALLEL SIMULATIONS (MULTIPLE SCENARIOS)  
  #              
  ###################################################

 
# To run multiple scenarios with varying parameters, use the #purrr:pmap" function


# create the different simulation scenarios with unique ids
sim_scenarios <- expand.grid(
  init_frequency = c(0.05, 0.1, 0.25, 0.5),
  # fecundity = c(1, 5, 10, 20),
  # n_loci = c(5,10,20,50,100),
  dispersal_prob = c(0.00001, 0.0001, 0.0005, 0.001)
) %>%
  mutate(
    scenario_id = row_number(),
    sim_id = paste0("freq", init_frequency, 
                    "_disp", dispersal_prob)
  )


# create a storgae folder for simulation results 
if (!dir.exists("results")) dir.create("results")


# progress bar to track simulation status
handlers(global = TRUE)
with_progress({
  sim_bar <- progressor(steps = nrow(sim_scenarios))
  
  
  # simulation
  sim <- sim_scenarios %>%
    pmap(function(init_frequency, dispersal_prob, scenario_id, sim_id) {
      sim_out <- run_model(
        n_patches = patches,
        pop_patches = pop_patches,
        n_per_patch = n_per_patch,
        n_loci = n_loci,
        init_frequency = init_frequency,
        fecundity = fecundity,
        carrying_capacity = carrying_capacity,
        prob_survival = prob_survival,
        decay = decay,
        lethal_effect = FALSE,
        complete_sterile = TRUE,
        sim_years = sim_years,
        overlapping = TRUE,
        adjacency_matrix = TRUE,
        dispersal_frac = dispersal_prob,
        decay
      )
      
      file_name <- paste0("results/sim_", sim_id, "_ID", scenario_id, ".rds")
      saveRDS(sim_out, file_name)
      
      sim_bar()
      file_name
    })
})




# # To run multiple scenarios with varying parameters, use the #purrr:pmap" function
# 
# 
# # create the different simulation scenarios with unique ids
# sim_scenarios <- expand.grid(
#   init_frequency = c(0.05, 0.1, 0.25, 0.5),
#   fecundity = c(1, 5, 10, 20),
#   n_loci = c(5,10,20,50,100),
#   dispersal_prob = c(0.00001, 0.0001, 0.0005, 0.001)
# ) %>%
#   mutate(
#     scenario_id = row_number(),
#     sim_id = paste0("freq", init_frequency, 
#                     "_fec", fecundity, 
#                     "_loci", n_loci, 
#                     "_disp", dispersal_prob)
#     )
# 
# 
# # create a storgae folder for simulation results 
# if (!dir.exists("results")) dir.create("results")
# 
# 
# # progress bar to track simulation status
# handlers(global = TRUE)
# with_progress({
#   sim_bar <- progressor(steps = nrow(sim_scenarios))
# 
#   
# # simulation
# sim <- sim_scenarios %>%
#   pmap(function(init_frequency, fecundity, n_loci, dispersal_prob, scenario_id, sim_id) {
#     sim_out <- run_model(
#       n_patches = patches,
#       pop_patches = pop_patches,
#       n_per_patch = n_per_patch,
#       n_loci = n_loci,
#       init_frequency = init_frequency,
#       fecundity = fecundity,
#       carrying_capacity = carrying_capacity,
#       prob_survival = prob_survival,
#       decay = decay,
#       lethal_effect = FALSE,
#       complete_sterile = TRUE,
#       sim_years = sim_years,
#       overlapping = TRUE,
#       adjacency_matrix = TRUE,
#       dispersal_frac = dispersal_prob,
#       decay
#     )
#     
#     file_name <- paste0("results/sim_", sim_id, "_ID", scenario_id, ".rds")
#     saveRDS(sim_out, file_name)
#     
#     sim_bar()
#     file_name
#     })
#   })
# 
# 
# 
