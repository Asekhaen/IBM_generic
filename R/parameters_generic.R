
# Parameters
growth_rate <- 1.5
carry_capacity <- 10000
prob_survival <- 0.76
dispersal_rate <- 0.103
patches <- 10                              # Number of patches
n_per_patch <- c(10000,0,0,0,0,0,0,0,0,0)            # Initial number of individuals per patch
growth_rate <- 2.5                            # Number of offspring per day per female mosquito
beta <- 100                               # the adult male population size at which the daily probability of mating is 0.5.
sim_days <-100                             # Number of simulation in days

# dispersal parameters
lambda <- 0.1
dispersal_frac <- 0.02


# Genetic (load) & drive parameters
n_loci <- 5
init_frequency = 0.25                    # initial frequency of deleterious recessives
decay <- 0.5  
