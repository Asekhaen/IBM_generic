
# Parameters
fecundity <- 5
density_dependence_factor <- 0.001
prob_survival <- 0.76
dispersal_rate <- 0.103
patches <- 25                               # Number of patches
n_per_patch <- c(1000,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0)             # Initial number of individuals per patch
growth_rate <- 2.5                          # Number of offspring per day per female mosquito
beta <- 100                                 # the adult male population size at which the daily probability of mating is 0.5.
sim_days <-20                               # Number of simulation in days

# dispersal parameters
lambda <- 0.1
dispersal_frac <- 0.02


# Genetic (load) & drive parameters
n_loci <- 20
init_frequency <- 0.25                   
decay <- 0.5  




# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- (runif(loci, max = genome.size))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}

l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)