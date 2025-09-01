
# Parameters
fecundity <- 5                             # Number of offspring per day per female mosquito
max_survival <- 1
prob_survival <- 0.7
decay_rate <- 2
dd_rate <- 0.0001
patches <- 20                              # Number of patches
n_per_patch <- c(40000,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,0)    # Initial number of individuals per patch
sim_years <- 50                            # Number of simulation in days
carrying_capacity = 40000                  # carrying capacity  

# dispersal parameters
lambda <- 0.1
dispersal_prob <- 0.0001

# Genetics: load/drive parameters
n_loci <- 20
init_frequency <- 0.25                   
decay <- 0.5  





reps <- 10
