
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


source("dependencies.R")

results <- run_model (
  n_patches = patches,
  pop_patches,
  n_per_patch = n_per_patch,
  n_loci = n_loci,
  init_frequency = init_frequency,
  fecundity = fecundity,
  carrying_capacity = carrying_capacity,
  prob_survival = prob_survival,
  dd_rate = dd_rate,
  decay = decay,
  lambda = lambda,
  lethal_effect = FALSE,
  complete_sterile = TRUE,
  sim_years = sim_years,
  adjacency_matrix = TRUE,
  dispersal_frac = dispersal_prob
)


if (!dir.exists("output")) dir.create("output")

# Save the results object inside the output folder
saveRDS(results, file = file.path("output", "results.rds"))


# # running replicates 
# 
# all_patch_stats <- list()
# all_allele_frequency <- list()
# 
# for (rep in 1:n_replicates) {
# results <- run_model (
#     n_patches = patches,
#     pop_patches,
#     n_per_patch = n_per_patch,
#     n_loci = n_loci,
#     init_frequency = init_frequency,
#     fecundity = fecundity,
#     carrying_capacity = carrying_capacity,
#     prob_survival = prob_survival,
#     dd_rate = dd_rate,
#     decay = decay,
#     lambda = lambda,
#     lethal_effect = FALSE,
#     complete_sterile = FALSE,
#     sim_years = sim_years,
#     adjacency_matrix = TRUE,
#     dispersal_frac = dispersal_prob
#   )
# 
# # Append to collectors
# all_patch_stats[[rep]] <- results$pop_stats
# all_allele_frequency [[rep]] <- results$allele_freq_per_locus
# }
# 
# 
# # append all replicates
# 
# # pop size 
# all_patch_stats_df <- bind_rows(
#   lapply(seq_along(all_patch_stats), function(rep) {
#     df <- all_patch_stats[[rep]]
#     df$replicate <- rep   # add replicate column
#     df
#   })
# )
# 
# #allele frequence per locus
# all_allele_frequency_df <- dplyr::bind_rows(
#   lapply(seq_along(all_allele_frequency), function(rep) {
#     df <- as.data.frame(all_allele_frequency[[rep]])
#     df$replicate <- rep
#     df
#   })
# )
# 
# saveRDS(all_patch_stats_df, file = "C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\100_rep_normal.rds")
# saveRDS(all_allele_frequency_df, file = "C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\100_normal_allele.rds")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# if (!dir.exists("output")) dir.create("output")
# saveRDS(results, file = "C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\results.rds")
# output1 <- (readRDS("C:\\Users\\22181916\\Documents\\Curtin-PhD\\R_and_IBM\\Generic_IBM_Proj\\IBM_generic\\output\\results.rds"))


  
  ###################################################
  #
  #   RUN PARALLEL SIMULATIONS (MULTIPLE SCENARIOS)  
  #              
  ###################################################

 
# To run multiple scenarios with varying parameters, use the #purrr:pmap" function



# # create the different simulation scenarios with unique ids
# sim_scenarios <- expand.grid(
#   #init_frequency = c(0.05, 0.1, 0.25),
#   lethal_effect = c(TRUE, FALSE),
#   complete_sterile = c(TRUE, FALSE),
#   dispersal_prob = c(0.00001, 0.0001, 0.001)
# 
# ) %>%
#   mutate(
#     scenario_name = paste0(#"freq", init_frequency,
#                            #"_fec", fecundity,
#                            #"_loci", n_loci,
#                            "_lethal", lethal_effect,
#                            "_sterle", complete_sterile,
#                            "_disp", dispersal_prob),
#     sim_id = row_number()
#   )
# 
# 
# sim_scenarios <- sim_scenarios[c(2,3,4,6,7,8,10,11,12),]
# 
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
#   # simulation
#   sim <- sim_scenarios %>%
#     pmap(function(dispersal_prob,
#                   scenario_name,
#                   lethal_effect,
#                   complete_sterile,
#                   sim_id) {
#       sim_out <- run_model(
#         n_patches = patches,
#         pop_patches,
#         n_per_patch = n_per_patch,
#         n_loci = n_loci,
#         init_frequency = init_frequency,
#         fecundity = fecundity,
#         carrying_capacity = carrying_capacity,
#         prob_survival = prob_survival,
#         dd_rate = dd_rate,
#         decay = decay,
#         lambda = lambda,
#         lethal_effect = TRUE,
#         complete_sterile = FALSE,
#         sim_years = sim_years,
#         overlapping = TRUE,
#         adjacency_matrix = TRUE,
#         dispersal_frac = dispersal_prob
#       )
# 
#       file_name <- paste0("results/sim_", sim_id, "_", scenario_name, ".rds")
#       saveRDS(sim_out, file_name)
# 
#       sim_bar()
#       file_name
#     })
# })
# 
# 






# # create the different simulation scenarios with unique ids
# sim_scenarios <- expand.grid(
#   init_frequency = c(0.05, 0.1, 0.25),
#   fecundity = c(1, 5, 10),
#   n_loci = c(10,50,100),
#   dispersal_prob = c(0.00001, 0.0001, 0.001)
# ) %>%
#   mutate(
#     scenario_name = paste0("freq", init_frequency,
#                            "_fec", fecundity,
#                            "_loci", n_loci,
#                            "_disp", dispersal_prob),
#     sim_id = row_number()
#   )

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
#   # simulation
#   sim <- sim_scenarios %>%
#     pmap(function(init_frequency,
#                   fecundity,
#                   n_loci,
#                   dispersal_prob,
#                   scenario_name,
#                   sim_id) {
#       sim_out <- run_model(
#         n_patches = patches,
#         pop_patches,
#         n_per_patch = n_per_patch,
#         n_loci = n_loci,
#         init_frequency = init_frequency,
#         fecundity = fecundity,
#         carrying_capacity = carrying_capacity,
#         prob_survival = prob_survival,
#         dd_rate = dd_rate,
#         decay = decay,
#         lambda = lambda,
#         lethal_effect = FALSE,
#         complete_sterile = TRUE,
#         sim_years = sim_years,
#         overlapping = TRUE,
#         adjacency_matrix = TRUE,
#         dispersal_frac = dispersal_prob
#       )
# 
#       file_name <- paste0("results/sim_", sim_id, "_", scenario_name, ".rds")
#       saveRDS(sim_out, file_name)
# 
#       sim_bar()
#       file_name
#     })
# })
