

# Density-dependent reproduction function (add a reference here)

bev_holt <- function(n_pop, fecundity, carrying_capacity) {
  return(fecundity / (1 + ((fecundity - 1) / carrying_capacity) * n_pop))
}



# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- sort((runif(loci, max = genome.size)))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}

# random selection of allele, with linkage

which_allele_fn <- function(n_offspring, n_loci, cov_matrix){
  epsilon <- MASS::mvrnorm(n_offspring, rep(0, n_loci), Sigma = cov_matrix)
  selection_prob <- plogis(epsilon)
  matrix(rbinom(n_offspring * n_loci, 1, selection_prob) == 1,
         nrow = n_offspring,
         ncol = n_loci)
}

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


# negative exponential dispersal kernel  

metapop <- function(coords, lambda, dispersal_frac) {
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
  
  return(dispersal_matrix)
}




# adjacency matrix 

step_stone <- function(n_patches, dispersal_frac) {
  
  matrix_landscape <- matrix(0, n_patches, n_patches)
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
  
  return(dispersal_matrix)
}

