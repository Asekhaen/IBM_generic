
# Parameters
fecundity <- 5                              # Number of offspring per day per female mosquito
prob_survival <- 0.7
dd_rate <- 0.0001
patches <- 7                               # Number of patches
n_per_patch <- c(10000,0,0,0,0,0,0)    # Initial number of individuals per patch
beta <- 100                           # the adult male population size at which the daily probability of mating is 0.5.
sim_years <- 5                         # Number of simulation in days
carrying_capacity = 10000             # carrying capacity  

# dispersal parameters
lambda <- 0.1
dispersal_prob <- 0.0005

# Genetics: load/drive parameters
n_loci <- 10
init_frequency <- 0.25                   
decay <- 0.5  

