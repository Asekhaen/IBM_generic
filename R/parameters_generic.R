
# Parameters
fecundity <- 5
density_dependence_factor <- 0.001
prob_survival <- 0.76
dispersal_rate <- 0.00103
patches <- 10                               # Number of patches
n_per_patch <- c(1000,0,0,0,0,0,0,0,0,0)    # Initial number of individuals per patch
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


# create coordinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


# create a dispersal matrix using the created function 
dispersal_matrix <- make_dispersal_matrix(coords = coords, 
                                          lambda = lambda, 
                                          dispersal_frac = dispersal_frac)

l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)
