

# population setup: Initialisation ####
ini_pop <- function(n_patches, n_per_patch, n_loci, init_frequency) {
  patches_pop <- list()
  for (i in 1:n_patches) {
    patches_pop[[i]] <- tibble(
      allele1 = matrix(rbinom(n = n_per_patch[i] * n_loci, 
                              size = 1, prob = init_frequency), ncol = n_loci), # 0 = wild-type, 1 = drive allele
      allele2 = matrix(rbinom(n = n_per_patch[i] * n_loci, 
                              size = 1, prob = init_frequency), ncol = n_loci),
      mate_allele1 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      mate_allele2 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      alive = TRUE
    )
    if (length(n_per_patch) != n_patches) warning("Initial patch population does not equal specified number of patches")
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
                   dd_rate,
                   overlapping,
                   sim_years) {
  #browser()
  #if(sim_years == 5) browser()
  updated_pop_patches <- list()
  updated_frac_homo <- list()
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    # reproduction  
    n.pop <- nrow(pop)
    
    if (n.pop > 0){
      
      
      # exp_fecundity <- fec_dd(n.pop, dd_rate, prob_survival)
      # act_fecundity <- rpois(n.pop, exp_fecundity)

      # exp_fecundity <- bev_holt(n.pop, dd_rate, fecundity)
      # act_fecundity <- rpois(n.pop, exp_fecundity)

      exp_fecundity <- bev_holt(n.pop, fecundity, carrying_capacity)
      act_fecundity <- rpois(n.pop, exp_fecundity)

      selected_mate_idx <- sample(n.pop, n.pop, replace = TRUE)
      selected_mate <- pop[selected_mate_idx,]
      pop$mate_allele1 <- selected_mate$allele1
      pop$mate_allele2 <- selected_mate$allele2
      
      
      if (complete_sterile) {
        # homozygous <- (homo_loci > 0)
        #homozygous <- rowSums((pop$allele1 + pop$allele2) == 2)        # homozygous loci for individual parent
        
        # homozygous <- rowSums(                                       # if either parent is homozygous
        #   ((pop$allele1 + pop$allele2) == 2) |
        #     ((pop$mate_allele1 + pop$mate_allele2) == 2)
        # )

        homozygous <- rowSums(
          ((pop$allele1 + pop$allele2) == 2) |
            ((pop$mate_allele1 + pop$mate_allele2) == 2)
        ) > 0
        
        sterile <- as.numeric(!homozygous)
        n_offspring <- act_fecundity * sterile
        
       # calculate the fraction that are homozygous for reproductive load (genetic load)
        n_homo <- sum(sterile == 0)
        n_individual <- length(sterile)
        frac_homo <- ifelse(n_individual > 0, n_homo / n_individual, 0)
        
      } else {
        n_offspring <- act_fecundity
        frac_homo <- 0
      }
      
    }   else {
      # If not, set offspring count to 0
      n_offspring <- rep(0, n.pop)
      frac_homo <- 0
      
    }
    
    
    # Offspring generation: Draw the actual number of offspring from a Poisson distribution
    
    # n_offspring <- rpois(n.pop, n_offspring)
    
    total_offspring <- sum(n_offspring)
    
    
    if (total_offspring > 0){  
      # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information
      
      ind_germline <- pop[rep(1:n.pop, n_offspring), c("allele1", "allele2")]
      mate_germline <- pop[rep(1:n.pop, n_offspring), c("mate_allele1", "mate_allele2")]
      

      # create covariance matrix betwee loci
      cov_matrix <- place_loci_mat(n_loci, genome.size = 1, var = 1, decay)
      
      
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
      
   
      # Update pop with offspring 
      pop <- offspring
      # pop <- bind_rows(pop, offspring)
     
  
    } else{
      pop <- tibble(
        allele1 = matrix(numeric(0), nrow = 0, ncol = n_loci),
        allele2 = matrix(numeric(0), nrow = 0, ncol = n_loci),
        mate_allele1 = matrix(numeric(0), nrow = 0, ncol = n_loci),
        mate_allele2 = matrix(numeric(0), nrow = 0, ncol = n_loci),
        alive = logical(0)
      )}
    # if turned on, this "if" statement simulates lethal effect of for individuals with 
    # homologous deleterious allele
    
    if (lethal_effect){
      homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
      any_homozygous <- rowSums(homozygous_lethal) > 0
      #pop <- filter(pop, !any_homozygous)
      pop <- pop[!any_homozygous,]
    }
    updated_pop_patches[[i]] <- pop
    updated_frac_homo[[i]] <- frac_homo
  }
  #return(updated_pop_patches)
  result <- list(
    updated_pop_patches = updated_pop_patches,
    updated_frac_homo = updated_frac_homo
  )
  
  return(result)
}


dispersal <- function(pop, n_patches, lambda, dispersal_frac, adjacency_matrix, check = TRUE) {
 #browser()
  
  # create a dispersal matrix using the created function 
  disp_matrix <- create_dispersal_matrix(n_patches, lambda, dispersal_frac, adjacency_matrix)
  
  patch_indices <- dispersed_pop <- vector(mode = "list", length = nrow(disp_matrix))
  
  # get new patch indices for all individuals
  for (i in seq_along(pop)) {
    patch <- pop[[i]]
    n_individuals<- nrow(patch)
    
    if (n_individuals == 0) next  # skip patch if no individuals
    
    dispersal_probs <- disp_matrix[i, ]   #dispersal probability for each patch based on the kernel matrix 
    
    # Desitnation patch for each individuals in this patch
    destination_patch <- sample(1:length(dispersal_probs), size = n_individuals, replace = TRUE, prob = dispersal_probs)

    # dataframe holding individuals and their destination patch
    patch_indices[[i]] <- tibble(ind_ID = seq_len(n_individuals),
                                 destination_patch = destination_patch)
  }
  
  #browser()
  # Move individuals to new patches
  for (i in seq_along(pop)) {
    patch <- pop[[i]]
    dispersers <- patch[patch_indices[[i]]$ind_ID, ]
    
    for (destination in seq_along(pop)) {
      ind_for_this_destination <- dispersers[patch_indices[[i]]$destination_patch == destination, ]
      dispersed_pop[[destination]] <- bind_rows(dispersed_pop[[destination]], ind_for_this_destination)
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
run_model <- function(n_patches,
                      pop_patches,
                      n_per_patch,
                      n_loci,
                      init_frequency,
                      fecundity,
                      carrying_capacity,
                      lethal_effect,
                      complete_sterile,
                      sim_years,
                      prob_survival,
                      dd_rate,
                      overlapping,
                      lambda, 
                      adjacency_matrix,
                      dispersal_frac,
                      decay
                      ){
  
  pop <- ini_pop(n_patches, n_per_patch, n_loci, init_frequency)
  
  patch_stats <- list()
  allele_frequency <- list()
  allele_freq_per_locus <- list()

  for (year in 1:sim_years) {
    #browser()
    cat("year", year, "Underway \n")
    
    # Growth with reproduction
    grown_pop <- growth(pop_patches = pop, 
                  n_loci,
                  carrying_capacity,
                  fecundity,
                  lethal_effect,
                  complete_sterile,
                  prob_survival,
                  dd_rate,
                  overlapping,
                  sim_years = year)
    
    pop <- grown_pop$updated_pop_patches
    frac_homo <- grown_pop$updated_frac_homo
    
    # Dispersal
    pop <- dispersal(pop, 
                     n_patches, 
                     lambda, 
                     dispersal_frac, 
                     adjacency_matrix, 
                     check = FALSE)

    if (nrow(pop[[n_patches]]) > carrying_capacity/2) {
      break
    }
    
    # Track annual population sizes per patch, occupancy rates, etc.

    patch_stats[[year]] <- tibble(
      year = year,
      patch = seq_along(pop),
      pop_size = sapply(pop, nrow),
      patch_occupied = sum(pop_size > establish_threshold),
      unoccupied = length(patch) - patch_occupied,
      prop_homo = unlist(frac_homo),
      occupancy_rate = patch_occupied/length(patch)
    )
  
    
    # Track proportion homozygous 

    # Track annual overall allele frequency per patch
    
  allele_frequency[[year]] <-  lapply(seq_along(pop), function(patch_id) {
    patch_pop <- pop[[patch_id]]
    loci_n  <- ncol(patch_pop$allele1)  
    n_ind   <- nrow(patch_pop$allele1)  
    total_allele_overall <- 2 * n_ind * loci_n
    
    overall <- tibble(
      year        = year,
      patch      = patch_id,
      total      = total_allele_overall,
      deleterious= sum(patch_pop$allele1 == 1) + sum(patch_pop$allele2 == 1),
      wild       = total_allele_overall - deleterious,
      freq       = ifelse(total_allele_overall == 0, 0, deleterious / total_allele_overall)
    )
  })
  
  
  # allele frequency per locus per patch
  
  allele_freq_per_locus[[year]] <- lapply(seq_along(pop), function(patch_id) {
    patch_pop <- pop[[patch_id]]
    loci_n    <- ncol(patch_pop$allele1)
    n_ind     <- nrow(patch_pop$allele1)
    total_per_locus <- 2 * n_ind        # per locus = 2 alleles per individual
    
    # deleterious alleles per locus
    deleterious_per_locus <- colSums(patch_pop$allele1 == 1) +
      colSums(patch_pop$allele2 == 1)
    
    # calculate frequencies per locus 
    freq_per_locus <- if (total_per_locus == 0) {
      rep(0, loci_n)
    } else {
      deleterious_per_locus / total_per_locus
    }
    
    tibble(
      year  = year,
      patch = patch_id,
      !!!setNames(as.list(freq_per_locus), paste0("locus", seq_len(loci_n)))
    )
  })
  
  #bind population dynamics outputs 
  patch_stats_df <- bind_rows(patch_stats)
  # calculate some stats
  patch_stats_df <- compute_stat(patch_stats_df,establish_threshold, carrying_capacity)
  # estimate spread rate or speed
  patch_stats_df <- invasion_speed(patch_stats_df)
  
  #bind combine allele frequency, select and add to population dynamics data
  allele_frequency_df <- bind_rows(allele_frequency)
  allele_freq <- allele_frequency_df |> select(freq)
  pop_stat <- bind_cols(patch_stats_df, allele_freq)
  
  #calculate fitness and genetic load
  

  #bind per locus allele frequency
  allele_freq_per_locus_df <- bind_rows(allele_freq_per_locus)
  
  }
  
  
  # Return the collected data
  results <- list(
    #final_pop = pop,
    pop_stats = pop_stat,
    allele_freq_per_locus = allele_freq_per_locus_df
  )
}

