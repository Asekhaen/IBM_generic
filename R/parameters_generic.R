
# Parameters
fecundity <- 5                             # Number of offspring per day per female mosquito
#max_survival <- 1
prob_survival <- 0.75
#decay_rate <- 2
dd_rate <- 0.004
patches <- 5                              # Number of patches
n_per_patch <- c(1000,0,0,0,0)  # Initial number of individuals per patch
sim_years <- 30                            # Number of simulation in days
carrying_capacity = 1000                  # carrying capacity  
establish_threshold <- 0.001 * carrying_capacity # 0.1% of carrying capacity 


# dispersal parameters
lambda <- 0.1
dispersal_prob <- 0.001

# Genetics: load/drive parameters
n_loci <- 10
init_frequency <- 0.25                   
decay <- 0.5  


# sensitivity analysis
n_samples <- 20
n_replicates <- 5