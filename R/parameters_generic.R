
# Parameters
fecundity <- 3                             # Number of offspring per day per female mosquito
#max_survival <- 1
prob_survival <- 0.75
#decay_rate <- 2
dd_rate <- 0.004
patches <- 2                              # Number of patches
carrying_capacity = (fecundity - 1)/dd_rate                  # carrying capacity  
n_per_patch <- c(carrying_capacity,0)             # Initial number of individuals per patch for two patch simulation
#n_per_patch <- c(carrying_capacity,0,0,0,0)             # Initial number of individuals per patch
sim_years <- 100                          # Number of simulation in days
establish_threshold <- round(0.005 * carrying_capacity) # 0.5% of carrying capacity 


# dispersal parameters
lambda <- 0.1
dispersal_prob <- 0.002

# Genetics: load/drive parameters
n_loci <- 1000                              # try varying 150 to 200
init_frequency <- 0.007                   
decay <- 0.5  

n_replicates <- 100
n_samples <- 20
