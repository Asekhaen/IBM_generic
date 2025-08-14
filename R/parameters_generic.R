
# Parameters
fecundity <- 5                              # Number of offspring per day per female mosquito
prob_survival <- 0.7
patches <- 7                               # Number of patches
n_per_patch <- c(10000,0,0,0,0,0,0)    # Initial number of individuals per patch
beta <- 100                           # the adult male population size at which the daily probability of mating is 0.5.
sim_days <-20                         # Number of simulation in days
carrying_capacity = 10000             # carrying capacity  

# dispersal parameters
lambda <- 0.1
dispersal_prob <- 0.005


# Genetic (load) & drive parameters
n_loci <- 10
init_frequency <- 0.25                   
decay <- 0.5  


# create coordinates for the patches/locations 
coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")


# create a dispersal matrix using the created function 
neg_exponet_model <- metapop (coords = coords, 
                              lambda = lambda, 
                              dispersal_frac = dispersal_prob)


adjacency_matrix <- step_stone(n_patches = patches, 
                               dispersal_frac = dispersal_prob)


l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)
