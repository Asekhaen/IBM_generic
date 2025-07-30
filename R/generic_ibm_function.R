library(tidyverse)
library(dplyr)
library(ggplot2)



# Beverton-Holt function
# soft density-dependence
# bev_holt <- function(N_pop, growth_rate, carry_capacity) {
#   return((growth_rate * N_pop) / (1 + (growth_rate - 1) * N_pop / carry_capacity))
# }

# hard density-dependence
bev_holt <- function(N_pop, growth_rate, carry_capacity) {
  return((growth_rate * N_pop) / (1 + N_pop / carry_capacity))
}




coords <- as.data.frame(100 * matrix(runif(patches * 2), ncol = 2))
colnames(coords) <- c("x","y")



ini_pop <- function(patches, n_per_patch, coord, n_loci, init_frequency) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      allele1 = matrix(rbinom(n = n_per_patch[i] * n_loci, size = 1, prob = init_frequency), ncol = n_loci), # 0 = wild-type, 1 = drive allele
      allele2 = matrix(rbinom(n = n_per_patch[i] * n_loci, size = 1, prob = init_frequency), ncol = n_loci),
      mate_allele1 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      mate_allele2 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      alive = TRUE
    )
    if (length(n_per_patch) != patches) warning("Initial patch population does not equal specified number of patches")
  }
  
  return(patches_pop)
}  


# pop <- ini_pop(patches,
#                n_per_patch,
#                n_loci = n_loci,
#                init_frequency = init_frequency)


# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- (runif(loci, max = genome.size))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}

l.cov.mat <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)




growth <- function(pop_patches, 
                   n_loci,
                   carry_capacity,
                   growth_rate,
                   lethal_effect,
                   complete_sterile,
                   sim_days,
                   loci_cov_matrix
                   ) {
  # if (sim_days == 15) browser()
  # browser()
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    n.pop <- nrow(pop)
    
    # growth 
    if (n.pop > 0){
      
      grown_pop <- bev_holt(n.pop, growth_rate = growth_rate, carry_capacity)
      avg_fecundity <- grown_pop/n.pop
      act_fecundity <- rpois(n.pop, avg_fecundity)
      
      selected_mate_idx <- sample(n.pop, n.pop, replace = TRUE)
      selected_mate <- pop[selected_mate_idx,]
      pop$mate_allele1 <- selected_mate$allele1
      pop$mate_allele2 <- selected_mate$allele2
      
      
      if (complete_sterile) {
        # homozygous <- (homo_loci > 0)
        homozygous <- rowSums((pop$allele1 + pop$allele2) == 2)        # homozygous loci for each female
        sterile <- as.numeric(!homozygous)
        n_offspring <- act_fecundity * sterile
      } else {
        n_offspring <- act_fecundity
      }
      
    }   else {
      # If not, set offspring count to 0
      n_offspring <- rep(0, n.pop)
    }
    
    
    # Offspring generation: Draw the actual number of offspring from a Poisson distribution
    
    total_offspring <- sum(n_offspring)
    
    
    if (total_offspring > 0){  
      # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information
      
      ind_germline <- pop[rep(1:n.pop, n_offspring), ] |> select(contains("allele"))
      mate_germline <- pop[rep(1:n.pop, n_offspring), ] |> select(contains("mate_allele"))
      
      
      # Genetic inheritance
      num_loci <- ncol(ind_germline$allele1)
      stopifnot(num_loci == n_loci)
      
      # # random selection of allele, with linkage 
      which_allele_fn <- function(n_offspring, num_loci, loci_cov_matrix){
        epsilon <- MASS::mvrnorm(n_offspring, rep(0, num_loci), Sigma = loci_cov_matrix)
        selection_prob <- plogis(epsilon)
        matrix(rbinom(n_offspring * num_loci, 1, selection_prob) == 1,
               nrow = n_offspring,
               ncol = num_loci)
      }
      
      
      which_allele_ind <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # female gametes
      which_allele_mate <- which_allele_fn(total_offspring, num_loci, loci_cov_matrix) # male gametes
      
      #  Determination of offspring features
      offspring <- tibble(
        allele1 = ifelse(which_allele_ind,
                         ind_germline$allele1,
                         ind_germline$allele2),
        allele2 = ifelse(which_allele_mate,
                         mate_germline$mate_allele1,
                         mate_germline$mate_allele2),
        mate_allele1 = matrix(NA, ncol = n_loci),
        mate_allele2 = matrix(NA, ncol = n_loci),
        alive = TRUE
      )
      
      # Update pop with offspring & fem population
      pop <- offspring
    }
    
    
    # Genetic load: lethal effect
    
    if (lethal_effect){
      homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
      any_homozygous <- rowSums(homozygous_lethal) > 0
      pop <- filter(pop, !any_homozygous)
    }
    
    
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
}


# 
# grown_pop <- growth(pop_patches = pop,
#                     fecundity,
#                     n_loci,
#                     daily_survival = prob_survival,
#                     beta,
#                     lethal_effect = FALSE,
#                     complete_sterile = TRUE,
#                     loci_cov_matrix = l.cov.mat,
#                     carry_capacity = carry_capacity)



# Dispersal: the default is the metapopulation network. Otherwise this can be switched to a stepping stone model


#### Metapopulation dispersal function ####
# first make dispersal matrix

make_dispersal_matrix <- function(coords, lambda, dispersal_frac) {
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


# create a dispersal matrix
dispersal_matrix <- make_dispersal_matrix(coords = coords, 
                                          lambda = lambda, 
                                          dispersal_frac = dispersal_frac)

# make meta_population dispersal function

meta_dispersal <- function(pop, dispersal_matrix, check = FALSE) {
  
  patch_indices <- dispersed_pop <- vector(mode = "list", length = nrow(dispersal_matrix))
  
  # get new patch indices for all individuals
  for (i in seq_along(pop)) {
    patch <- pop[[i]]
    n_patch <- nrow(patch)
    
    if (n_patch == 0) next  # skip patch if no individuals
    
    dispersal_probs <- dispersal_matrix[i, ]
    
    # Sample new patch indices for all individuals in this patch
    new_pop_indices <- sample(1:length(dispersal_probs), size = n_patch, replace = TRUE, prob = dispersal_probs)
    
    # Store original row indices and their assigned new patch
    patch_indices[[i]] <- tibble(ind_within_pop_indices = seq_len(n_patch),
                                 new_pop_indices = new_pop_indices)
  }
  
  # Move individuals to new patches
  for (i in seq_along(pop)) {
    patch <- pop[[i]]
    dispersers <- patch[patch_indices[[i]]$ind_within_pop_indices, ]
    
    for (jj in 1:length(pop)) {
      ind_jj <- dispersers[patch_indices[[i]]$new_pop_indices == jj, ]
      dispersed_pop[[jj]] <- bind_rows(dispersed_pop[[jj]], ind_jj)
    }
  }
  
  if (check){
    n_pop <- sum(sapply(pop, nrow))
    n_disp <- sum(sapply(dispersed_pop, nrow))
    cat(n_pop, " ", n_disp, "\n")
  }
  
  return(dispersed_pop)
}


#### Stepping stone dispersal (one-dimensional discrete space) ####


ss_dispersal <- function(pop_patches, dispersal_frac) {
  n_patches <- length(pop_patches)
  dispersed_pop <- vector("list", n_patches)
  
  # # Initialize empty population for each patch
  # for (i in seq_len(n_patches)) {
  #   dispersed_pop[[i]] <- pop_patches[[i]][0, ]
  # }
  # 
  for (i in seq_len(n_patches)) {
    current_patch <- pop_patches[[i]]
    
    if (nrow(current_patch) == 0) next   # skip unoccupied patch
    
    # Determine which adults disperse
    ready_to_disperse <- rbinom(nrow(current_patch), 1, dispersal_frac)
    dispersers <- current_patch[ready_to_disperse == 1, ]
    non_dispersers <- current_patch[ready_to_disperse == 0, ]
    
    # Disperse adults to i-1 or i+1
    if (nrow(dispersers) > 0) {
      directions <- sample(c(-1, 1), nrow(dispersers), replace = TRUE)
      target_patch <- i + directions
      target_patch <- pmin(pmax(target_patch, 1), n_patches)  # keep within boundaries i.e. patch 1 and patch "n"
      
      for (j in seq_along(target_patch)) {
        dispersed_pop[[target_patch[j]]] <- bind_rows(
          dispersed_pop[[target_patch[j]]],
          dispersers[j, ]
        )
      }
    }
    
    # Add stayers (non-dispersing adults and non-adults) to current patch
    dispersed_pop[[i]] <- bind_rows(
      dispersed_pop[[i]], 
      non_dispersers
      )
  }
  
  return(dispersed_pop)
}



#### Simulation function 
simulation <- function(patches,
                       pop_patches,
                       n_per_patch,
                       n_loci,
                       init_frequency,
                       carry_capacity,
                       growth_rate,
                       lethal_effect = FALSE,
                       complete_sterile = TRUE,
                       sim_days = sim_days,
                       stepping_stone_model = FALSE,
                       loci_cov_matrix = l.cov.mat,
                       decay = decay,
                       dispersal_matrix = dispersal_matrix
                       ){
  
  pop <- ini_pop(patches,
                 n_per_patch,
                 n_loci = n_loci,
                 init_frequency = init_frequency)
  
  l.cov.mat <- place_loci_mat(n_loci, 
                              genome.size = 1, 
                              var = 1, decay)
  
  patch_sizes <- list()
  allele_frequency <- list()
  # spread_rate <- list()
  # gen_time <- list()
  
  
  for (day in 1:sim_days) {
    #if (day == 45) browser()
    cat("Day", day, "Underway \n")
    
    # Growth with reproduction
    pop <- growth(pop_patches = pop, 
                  n_loci,
                  carry_capacity,
                  growth_rate,
                  lethal_effect = FALSE,
                  complete_sterile = TRUE,
                  sim_days,
                  loci_cov_matrix = l.cov.mat)
    
    # Dispersal
    if(stepping_stone_model) {
      pop <- ss_dispersal(pop, dispersal_frac)
    } else {
      pop <- meta_dispersal(pop, dispersal_matrix, check = FALSE)
    }
    
    # 
    # # Track the average generation time: from egg to first oviposition
    # 
    # all_females[[day]] <- do.call(rbind, lapply(pop, function(pop) {
    #   filter(pop, sex == 1, !is.na(first_ovip_day), !is.na(birth))
    # }))
    # 
    # all_females$time_to_first_ovip <- all_females$first_ovip_day - all_females$birth
    # gen_time <- mean(all_females$time_to_first_ovip)
    
    
    # Track daily population sizes per patch
    
    patch_sizes[[day]] <- do.call(rbind, lapply(seq_along(pop), function(patch_id) {
      data.frame(
        day = day,
        patch = patch_id,
        pop_size = nrow(pop[[patch_id]])
      )
    }))
    patch_sizes_df <- do.call(rbind, patch_sizes)
    
    
    # Track daily allele frequency per patch
    
    allele_frequency[[day]] <- do.call(rbind, lapply(seq_along(pop), function(patch_id) {
      patch_pop <- pop[[patch_id]]
      if (nrow(patch_pop) > 0) {
        deleterious <- sum(patch_pop$allele1 == 1) + sum(patch_pop$allele2 == 1)
        total <- 2 * nrow(patch_pop) * ncol(patch_pop$allele1)
        wild_type <- total - deleterious
        freq <- deleterious / total
      } else {
        deleterious <- 0
        total <- 0
        wild_type <- 0
        freq <- 0
      }
      data.frame(
        patch = patch_id,
        wild = wild_type,
        lethal = deleterious,
        total = total,
        freq = freq,
        day = day
      )
    }
    ))
    allele_frequency_df <- do.call(rbind, allele_frequency)
  }
  
  # track spread or invasion rate
  
  # Return the collected data
  list(
    patch_sizes = patch_sizes_df,
    allele_frequency = allele_frequency_df,
    # spread_rate <-
    # gen_time <- gen_time,
    final_pop = pop
  )
}


