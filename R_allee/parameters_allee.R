
# Parameters
fecundity <- 3                             # Number of offspring per day per female mosquito
patches <- 6                              # Number of patches
carrying_capacity = 1000
half_K <- 0.75 * carrying_capacity
n_per_patch <- c(carrying_capacity, 5, 10, 20, 50, 100)             # Initial number of individuals per patch
sim_years <- 100                          # Number of simulation in days
colonise_threshold <- round(0.005 * carrying_capacity) # 0.5% of carrying capacity 


# dispersal parameters
lambda <- 0.1
dispersal_frac <- 0.0025

# Genetics: load/drive parameters
n_loci <- 1000                              # try varying 150 to 200
n_load <- 0.25
init_frequency = calc_q(n_load, n_loci)                                     # for analysis use 0.01, 0.025, 0.05, 0.1
decay <- 0.5  
n_replicates <- 25
n_samples <- 100





