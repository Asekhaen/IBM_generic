
# Parameters
fecundity <- 3                                            # Number of offspring per day per female mosquito
patches <- 25                                             # Number of patches
carrying_capacity = 1000                                  # carrying capacity
half_K <- 0.75 * carrying_capacity
n_per_patch <- create_n_per_patch(patches, 
                                  carrying_capacity)      # Initial number of individuals per patch
sim_years <- 200                                           # Number of generations 350 to 500 generations 
establish_threshold <- round(0.005 * carrying_capacity)   # 0.5% of carrying capacity 

# dispersal parameters
lambda <- 0.5                                             # dispersal decay parameter 
dispersal_frac <- 0.0025                                  # for analysis we used 0.001, 0.0025, 0.005, 0.01) 

# Genetics: load parameters
n_loci <- 1000                                             # for analysis use 1, 10, 100, 1000) 
init_frequency <- 0.01                                     # for analysis use 0.01, 0.025, 0.05, 0.1
# decay <- 0.5  

n_replicates <- 100                                        # make this 1000 replicates for the final (analysed) data set to capture more stochasticity
n_samples <- 100

