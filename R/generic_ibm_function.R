

# population setup: Initialisation ####
ini_pop <- function(patches, n_per_patch, n_loci, init_frequency) {
  patches_pop <- list()
  
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      allele1 = matrix(rbinom(n = n_per_patch[i] * n_loci, 
                              size = 1, prob = init_frequency), ncol = n_loci), # 0 = wild-type, 1 = drive allele
      allele2 = matrix(rbinom(n = n_per_patch[i] * n_loci, 
                              size = 1, prob = init_frequency), ncol = n_loci),
      mate_allele1 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      mate_allele2 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      alive = TRUE
    )
    if (length(n_per_patch) != patches) warning("Initial patch population does not equal specified number of patches")
  }
  
  return(patches_pop)
}  


# growth function ####
# captures density-dependent reproduction using Beverton-Holt model, genetic 
# (with linkage) inheritance, and whether the population is overlapping or 
# non-overlapping
# 

growth <- function(pop_patches, 
                   n_loci,
                   carrying_capacity,
                   fecundity,
                   lethal_effect,
                   complete_sterile,
                   prob_survival,
                   overlapping,
                   cov_matrix,
                   sim_days) {
  # if(sim_days == 5) browser()
  updated_pop_patches <- list()
  
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    n.pop <- nrow(pop)
    
    # growth 
    if (n.pop > 0){
      
      exp_fecundity <- bev_holt(n.pop, fecundity, carrying_capacity)
      act_fecundity <- rpois(n.pop, exp_fecundity)
      
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
      
      ind_germline <- pop[rep(1:n.pop, n_offspring), c("allele1", "allele2")]
      mate_germline <- pop[rep(1:n.pop, n_offspring), c("mate_allele1", "mate_allele2")]
      

      # Genetic inheritance
      
      which_allele_ind <- which_allele_fn(total_offspring, n_loci, cov_matrix) # female gametes
      which_allele_mate <- which_allele_fn(total_offspring, n_loci, cov_matrix) # male gametes
      
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
      
      if (overlapping) {
        pop <- bind_rows(pop, offspring)
        pop <- pop[rbinom(nrow(pop), 1, prob_survival) == 1, ]    # survival 
      } else{
        pop <- offspring
        pop <- pop[rbinom(nrow(pop), 1, prob_survival) == 1, ]    # survival 
      }
    }
    
    # if turned on, this "if" statement simulates lethal effect of for individuals with 
    # homologous deleterious allele
    
    if (lethal_effect){
      homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
      any_homozygous <- rowSums(homozygous_lethal) > 0
      #pop <- filter(pop, !any_homozygous)
      pop <- pop[pop[!any_homozygous],]
    }
    updated_pop_patches[[i]] <- pop
  }
  return(updated_pop_patches)
}

# dispersal: uses a negative exponential dispersal kernel (for spatial metapopulation) or  
# adjacency nearest neighbour (for one dimensional space stepping stone model)

# negative exponential kernel ####

dispersal <- function(pop, dispersal_type, check = FALSE) {
  patch_indices <- dispersed_pop <- vector(mode = "list", length = nrow(dispersal_type))
  
  # get new patch indices for all individuals
  for (i in seq_along(pop)) {
    patch <- pop[[i]]
    n_patch <- nrow(patch)
    
    if (n_patch == 0) next  # skip patch if no individuals
    
    dispersal_probs <- dispersal_type[i, ]
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

#### Simulation function 
run_model <- function(patches,
                      pop_patches,
                      n_per_patch,
                      n_loci,
                      init_frequency,
                      fecundity,
                      carrying_capacity,
                      lethal_effect,
                      complete_sterile,
                      sim_days,
                      prob_survival,
                      overlapping,
                      cov_matrix,
                      decay,
                      dispersal_type
                      ){
  
  pop <- ini_pop(patches,
                 n_per_patch,
                 n_loci,
                 init_frequency)
  

  patch_sizes <- list()
  allele_frequency <- list()
  # spread_rate <- list()
  # gen_time <- list()
  
  for (day in 1:sim_days) {
    #if (day == 5) browser()
    cat("Day", day, "Underway \n")
    
    # Growth with reproduction
    pop <- growth(pop_patches = pop, 
                  n_loci,
                  carrying_capacity,
                  fecundity,
                  lethal_effect,
                  complete_sterile,
                  prob_survival,
                  overlapping,
                  cov_matrix,
                  sim_days = day)
    
    # Dispersal
  
      pop <- meta_dispersal(pop, dispersal_type, check = FALSE)

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


