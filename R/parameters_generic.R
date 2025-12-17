
# Parameters
fecundity <- 3                             # Number of offspring per day per female mosquito
#max_survival <- 1
#prob_survival <- 0.75
#decay_rate <- 2
#dd_rate <- 0.004
patches <- 10                              # Number of patches
carrying_capacity = 1000                  # carrying capacity
half_K <- 0.5 * carrying_capacity
n_per_patch <- c(carrying_capacity,0,0,0,0,0,0,0,0,0)             # Initial number of individuals per patch for two patch simulation
#n_per_patch <- c(carrying_capacity,0,0,0,0)             # Initial number of individuals per patch
sim_years <- 5                          # Number of simulation in days
establish_threshold <- round(0.01 * carrying_capacity) # 0.5% of carrying capacity 

# dispersal parameters
lambda <- 3                   # dispersal decay parameter 
dispersal_frac <- 0.002

# Genetics: load/drive parameters
n_loci <- 175                 # try varying 150 to 200
init_frequency <- 0.035                   
decay <- 0.5  

n_replicates <- 5
n_samples <- 1000
