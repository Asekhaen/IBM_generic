# function to load packages 

load_libraries <- function(pack) {
  for (p in pack) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
}


# Density-dependent reproduction 

bev_holt <- function(n_pop, fecundity, carrying_capacity) {
  return(fecundity / (1 + (fecundity - 1) / carrying_capacity * n_pop))
}
 
# bev_holt <- function(n_pop, dd_rate, fecundity){
#   expected_offspring <- fecundity/(1 + dd_rate*n_pop)
#   return(expected_offspring)
# }

# Density dependent fecundity function
# fec_dd <- function(n_pop, dd_rate, prob_survival, c = 0.06){
#   fecundity <- (1-prob_survival)/c
#   expected_offspring <- fecundity*exp(-dd_rate*n_pop)
#   return(expected_offspring)
# }


# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  set.seed(11) # note that I added this to prevent the loci position from each simulation
  loci_positions <- sort((runif(loci, max = genome.size)))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}

# random selection of allele, with linkage
which_allele_fn <- function(n_offspring, n_loci, cov_matrix){
  epsilon <- MASS::mvrnorm(n_offspring, rep(0, n_loci), Sigma = cov_matrix)
  selection_prob <- epsilon > 0
  selection_prob[] <- as.numeric(selection_prob)
  matrix(rbinom(n_offspring * n_loci, 1, selection_prob) == 1,
         nrow = n_offspring,
         ncol = n_loci)
}


# which_allele_fn <- function(n_offspring, n_loci, cov_matrix){
#   epsilon <- MASS::mvrnorm(n_offspring, rep(0, n_loci), Sigma = cov_matrix)
#   selection_prob <- plogis(epsilon)
#   matrix(rbinom(n_offspring * n_loci, 1, selection_prob) == 1,
#          nrow = n_offspring,
#          ncol = n_loci)
# }

# alternative  function for computational speed
# which_allele_fn <- function(n_ind, n_loci, loci_cov_matrix) {
#   epsilon <- MASS::mvrnorm(n = n_ind,
#                            mu = rep(0, n_loci),
#                            Sigma = loci_cov_matrix)
#   
#   # alternatively, pass in 'L_loci_cov_matrix', which is computed earlier as:
#   #   L_loci_cov_matrix <- chol(loci_cov_matrix)
#   # then inside this function do:
#   #   z <- matrix(rnorm(n_ind * n_loci), n_ind, n_loci)
#   #   epsilon <- z %*% L
#   
#   selection_prob <- 1 / (1 + exp(-epsilon))
#   u <- matrix(runif(n_ind * n_loci),
#               n_ind, n_loci)
#   u < selection_prob
# }


# random selection of allele, without linkage
which_allele <- function(offspring, n_loci){
    matrix(rbinom(n_offspring * n_loci, 1, 0.5) == 1,
         nrow = n_offspring,
         ncol = n_loci)
}


# create random coordinates for the patches

create_coordinates <- function(x) {
  coords <- matrix(runif(x * 2), ncol = 2)
  return(coords)
}

# dispersal: uses a either a negative exponential dispersal kernel (for spatial
# metapopulation) or nearest neighbour/adjacency matrix (one dimensional space
# stepping stone model). This function  creates a dispersal kernel: either adjacency
# matrix or a negative exponential decay matrix based on a true of false statement
# (see run_file.R)

create_dispersal_matrix <- function(patches, lambda, dispersal_frac, adjacency_matrix){

  if(adjacency_matrix){
    matrix_landscape <- matrix(0, patches, patches)
    adjacency <- abs(row(matrix_landscape) - col(matrix_landscape)) == 1
    adjacency[] <- as.numeric(adjacency)

    # make these rows sum to 1 to get probability of moving to other patch
    # *if* they left. This dispersal matrix gives the probability of the vector
    # vector moving between patches
    rel_dispersal_matrix <- sweep(adjacency, 1,
                                  rowSums(adjacency), FUN = "/")

    # normalise these to have the overall probability of dispersing to that patch,
    # and add back the probability of remaining
    dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
      (1 - dispersal_frac) * diag(nrow(adjacency))

  } else {
    coords <- create_coordinates(patches)
    # dispersal matrix
    dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
    #exponential dispersal kernel
    dispersal_kernel <- exp(-lambda * dist_matrix)
    # set the diagonal elements to 0 to prevent self-dispersal
    diag(dispersal_kernel) <- 0

    # make these rows sum to 1 to get probability of moving to other patch
    # *if* they left. This dispersal matrix gives the probability of the vector
    # vector moving between patches
    rel_dispersal_matrix <- sweep(dispersal_kernel, 1,
                                  rowSums(dispersal_kernel), FUN = "/")

    # normalise these to have the overall probability of dispersing to that patch,
    # and add back the probability of remaining
    dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
      (1 - dispersal_frac) * diag(nrow(dispersal_kernel))


    # to ensure tat the probability of movement between patches aligns with the
    # number of individuals per patch when comparing the plot with the patch population statistics

  }
  return(dispersal_matrix)
}



# create_dispersal_matrix <- function(n_patches, lambda, dispersal_frac, adjacency_matrix){
#   
#   if(adjacency_matrix){
#     matrix_landscape <- matrix(0, n_patches, n_patches)
#     adjacency <- abs(row(matrix_landscape) - col(matrix_landscape)) == 1
#     adjacency[] <- as.numeric(adjacency)
#     
#     # make these rows sum to 1 to get probability of moving to other patch
#     # *if* they left. This dispersal matrix gives the probability of the vector
#     # vector moving between patches
#     rel_dispersal_matrix <- sweep(adjacency, 1,
#                                   rowSums(adjacency), FUN = "/")
#     
#     # normalise these to have the overall probability of dispersing to that patch,
#     # and add back the probability of remaining
#     dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
#       (1 - dispersal_frac) * diag(nrow(adjacency))
#   } else {
#     # create coordinates for the patches/locations
#     coords <- as.data.frame(100 * matrix(runif(n_patches * 2), ncol = 2))
#     colnames(coords) <- c("x","y")
#     
#     # dispersal matrix
#     dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
#     #exponential dispersal kernel
#     dispersal_kernel <- exp(-lambda * dist_matrix)
#     # set the diagonal elements to 0 to prevent self-dispersal
#     diag(dispersal_kernel) <- 0
#     # make these rows sum to 1 to get probability of moving to other patch
#     # *if* they left. This dispersal matrix gives the probability of the vector
#     # vector moving between patches
#     rel_dispersal_matrix <- sweep(dispersal_kernel, 1,
#                                   rowSums(dispersal_kernel), FUN = "/")
#     
#     # normalise these to have the overall probability of dispersing to that patch,
#     # and add back the probability of remaining
#     dispersal_matrix <- dispersal_frac * rel_dispersal_matrix +
#       (1 - dispersal_frac) * diag(nrow(dispersal_kernel))
#   }
#   return(dispersal_matrix)
# }



# function for invasion speed 
invasion_speed <- function(data) {
  occupied = data$patch_occupied
  time = data$year
  # fit model
  model <- lm(occupied ~ time)
  data$speed <- coef(model)[2]
  return(data)
}