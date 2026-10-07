#===============================================================
# GENERIC IBM MODEL (for invasive species population dynamics) 
#===============================================================

# =================================
# POPULATION INITIALISATION 
# =================================

ini_pop <- function(patches, n_per_patch, n_loci, init_frequency) {
  patches_pop <- list()
  for (i in 1:patches) {
    patches_pop[[i]] <- tibble(
      allele1 = matrix(rbinom(n = n_per_patch[i] * n_loci, 
                              size = 1, prob = init_frequency), ncol = n_loci), # 0 = wild-type, 1 = drive or deleterious allele
      allele2 = matrix(rbinom(n = n_per_patch[i] * n_loci, 
                              size = 1, prob = init_frequency), ncol = n_loci),
      mate_allele1 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      mate_allele2 = matrix(NA, nrow = n_per_patch[i], ncol = n_loci),
      alive = TRUE
    )
    if (length(n_per_patch) != patches) warning(
      "Initial patch population does not equal specified number of patches")
  }
  return(patches_pop)
}  

pop <- ini_pop(patches, n_per_patch, n_loci, init_frequency)

# ===============================================
# REPRODUCTION AND GENETIC INHERITANCE  
# ===============================================

# captures density-dependent reproduction using Beverton-Holt model, and genetics 
# (with linkage) inheritance 

growth <- function(pop_patches, 
                   n_loci,
                   carrying_capacity,
                   fecundity,
                   lethal_effect,
                   complete_sterile,
                   linkage,
                   sim_years) {
  # browser()
  # if(sim_years == 1) browser()
  updated_pop_patches <- list()
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    # reproduction  
    n.pop <- nrow(pop)
    
    if (n.pop >= 2){
      
      
      # # exp_fecundity <- fec_dd(n.pop, dd_rate, prob_survival)
      # # act_fecundity <- rpois(n.pop, exp_fecundity)

      exp_fecundity <- bev_holt(n.pop, fecundity, carrying_capacity)
      act_fecundity <- rpois(n.pop, exp_fecundity)
      
      # obligate outcrossing: each individual draws one mate uniformly from the other
      # n.pop - 1 individuals (random offset 1..n.pop-1 around the index, so self is never chosen)
      selected_mate_idx <- ((seq_len(n.pop) - 1 +
                               sample.int(n.pop - 1, n.pop, replace = TRUE)) %% n.pop) + 1
      
      # sterile = homozygous for the deleterious allele at >= 1 locus
      sterile <- rowSums((pop$allele1 + pop$allele2) == 2) > 0
      
      # sterile effect of homozygous deleterious alleles
      if (complete_sterile) {
        homo_del <- as.numeric(!(sterile | sterile[selected_mate_idx]))
        n_offspring <- act_fecundity * homo_del
      }  else {
        n_offspring <- act_fecundity
      }
      
    }   else {
      # 0 or 1 individual: no mate available, so no reproduction (obligate outcrossing)
      n_offspring <- rep(0, n.pop)
    }
    
     total_offspring <- sum(n_offspring)
     
    
 # GENETIC INHERITANCE
 # function with option to choose between inheritance type (linkage and without linkage)

     if (total_offspring > 0 && n.pop >= 2){   # mates need at least 2 individuals
       
       # attach each individual's mate genotype (mates drawn above)
       selected_mate <- pop[selected_mate_idx,]
       pop$mate_allele1 <- selected_mate$allele1
       pop$mate_allele2 <- selected_mate$allele2
       
       # Replicate the parents features `n_offspring` times for each offspring, collect only genetic information
       ind_germline <- pop[rep(1:n.pop, n_offspring), c("allele1", "allele2")]
       mate_germline <- pop[rep(1:n.pop, n_offspring), c("mate_allele1", "mate_allele2")]


       if (linkage) {

         # create covariance matrix between loci
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

       } else {
         # Generate random allele choices for each locus of each offspring
         which_allele_ind  <- matrix(rbinom(total_offspring * n_loci, 1, 0.5),
                                     nrow = total_offspring, ncol = n_loci)
         which_allele_mate <- matrix(rbinom(total_offspring * n_loci, 1, 0.5),
                                     nrow = total_offspring, ncol = n_loci)

         # Dataframe for Offspring
         offspring <- tibble(
           allele1 = ifelse(which_allele_ind,
                            ind_germline$allele1,
                            ind_germline$allele2),
           allele2 = ifelse(which_allele_mate,
                            mate_germline$mate_allele1,
                            mate_germline$mate_allele2),
           mate_allele1 = matrix(NA, ncol = n_loci),
           mate_allele2 = matrix(NA, ncol = n_loci),
           alive          = TRUE
         )
       }
       
       
       
       if (lethal_effect){
         homozygous_lethal <- (offspring$allele1 == 1) & (offspring$allele2 == 1)
         any_homozygous <- rowSums(homozygous_lethal) > 0
         # any_homozygous_del <- as.numeric(!any_homozygous)
         # n_homo_del <- sum(any_homozygous_del == 0)
         # n_individual <- length(any_homozygous_del)
         offspring <- offspring[!any_homozygous,]
    
       }
       
       # Update pop with offspring
       pop <- offspring

     } else{
       pop <- tibble(
         allele1 = matrix(numeric(0), nrow = 0, ncol = n_loci),
         allele2 = matrix(numeric(0), nrow = 0, ncol = n_loci),
         mate_allele1 = matrix(numeric(0), nrow = 0, ncol = n_loci),
         mate_allele2 = matrix(numeric(0), nrow = 0, ncol = n_loci),
         alive = logical(0)
       )}

    updated_pop_patches[[i]] <- pop
  }
  #return(updated_pop_patches)
  result <- list(
    updated_pop_patches = updated_pop_patches
  )
  
  return(result)
}


# ===========================
#  DISPERSAL
# ===========================


dispersal <- function(pop, patches, lambda, dispersal_frac, adjacency_matrix, check = TRUE) {
  disp_matrix <- create_dispersal_matrix(patches, lambda, dispersal_frac, adjacency_matrix)
  # create empty list of patches to hold dispersed pop also the column structure
  dispersed_pop <- vector("list", patches)
  for (i in seq_len(patches)) {
    dispersed_pop[[i]] <- pop[[i]][0, ]
  }
  # flows[i, j] = number of individuals moving from patch i to patch j (diagonal = stayers)
  flows <- matrix(0L, patches, patches)
  
  # extract individual from each patch and skip if there are no individual in the patch
  for (i in seq_along(pop)) {
    patch <- pop[[i]]
    n <- nrow(patch)
    if (n == 0) next
    
    #Extract dispersal probabilities from the adjacency or exponential dispersal krnel
    probs <- disp_matrix[i, ]
    
    # safety checks
    stopifnot(
      length(probs) == patches,
      all(probs >= 0),
      abs(sum(probs) - 1) < 1e-8
    )
    #sample the desitination for each individual
    destinations <- sample(
      seq_len(patches),
      size = n,
      replace = TRUE,
      prob = probs
    )
    flows[i, ] <- tabulate(destinations, nbins = patches)
    
    # bind individuals that moved to their destinbtion patch
    for (j in seq_len(patches)) {
      dispersed_pop[[j]] <- bind_rows(
        dispersed_pop[[j]],
        patch[destinations == j, ]
      )
    }
  }
  
  if (check) {
    cat(
      "Before:", sum(sapply(pop, nrow)),
      "After:",  sum(sapply(dispersed_pop, nrow)), "\n"
    )
  }
  
  list(pop = dispersed_pop, flows = flows)
  }

#================================ 
# SIMULATION FUNCTION  
#================================ 

run_model <- function(patches,
                      pop_patches,
                      n_per_patch,
                      n_loci,
                      init_frequency,
                      fecundity,
                      carrying_capacity,
                      lethal_effect,
                      complete_sterile,
                      linkage,
                      sim_years,
                      lambda, 
                      adjacency_matrix,
                      dispersal_frac,
                      decay
                      ){
  
  pop <- ini_pop(patches, n_per_patch, n_loci, init_frequency)
  
  patch_stats <- list()
  genetic_data <- list()
  time_halfK <- rep(NA, patches)
  dispersal_kernel <- matrix()
  
  # propagule supply / establishment tracking
  time_colonised        <- rep(NA_integer_, patches) # first year N >= colonisation_threshold (arrival time)
  cum_emigrants         <- rep(0, patches)           # cumulative individuals that left each patch
  cum_immigrants        <- rep(0, patches)           # cumulative individuals that arrived in each patch
  immigrants_to_arrival <- rep(NA_real_, patches)    # cumulative immigrants up to and including arrival year
  
  for (year in 1:sim_years) {
    #browser()
    cat("year", year, "Underway \n")

   # to track previous growth
    prev_pop_size <- sapply(pop, nrow)
    
    # Growth with reproduction
    grown_pop <- growth(pop, 
                  n_loci,
                  carrying_capacity,
                  fecundity,
                  lethal_effect,
                  complete_sterile,
                  linkage,
                  sim_years = year)
    
    pop <- grown_pop$updated_pop_patches

    # S_t: surviving offspring per patch before dispersal (after lethal removal)
    post_repro_size <- sapply(pop, nrow)
    
    #if (nrow(pop[[patches]]) > carrying_capacity/2) {
    #  break
    #}

    # Track genetic stats 
    genetic_data[[year]] <- lapply(seq_along(pop), function(patch_id) {
      #browser()
      patch_pop <- pop[[patch_id]]
      allele1   <- patch_pop$allele1   
      allele2   <- patch_pop$allele2   
      n_ind   <- nrow(allele1)
      n_loci  <- ncol(allele1)
      
      genotype_sum <- allele1 + allele2    
      AA_count <- colSums(genotype_sum == 0)   # homozygous wild/normal
      Aa_count <- colSums(genotype_sum == 1)   # heterozygous recessive
      aa_count <- colSums(genotype_sum == 2)   # homozygous deleterious (the proportion of "aa" can be use as measurement for genetic load)
      
      
      total_alleles <- 2 * n_ind
      deleterious_count <- colSums(allele1 == 1) + colSums(allele2 == 1)
      wild_count <- total_alleles - deleterious_count
      freq_a <- ifelse(deleterious_count > 0, deleterious_count / total_alleles, 0)
      freq_A <- ifelse(wild_count > 0, wild_count / total_alleles, 0)
      
      # genetic load (Eq. 4a): proportion of individuals homozygous deleterious at >= 1 locus
      realised_load <- if (n_ind > 0) mean(rowSums(genotype_sum == 2) > 0) else NA_real_
      # expected load under HWE and independent loci (Eq. 4), for comparison
      hwe_load <- if (n_ind > 0) 1 - prod(1 - freq_a^2) else NA_real_
      
      tibble::tibble(
        patch = patch_id,
        year  = year,
        locus = 1:n_loci,
        AA = AA_count,
        Aa = Aa_count,
        aa = aa_count,
        freq_A = freq_A,
        freq_a = freq_a,
        n_ind = n_ind,
        realised_load = realised_load,
        hwe_load = hwe_load
      )
    })
    
    #Dispersal
    disp <- dispersal(pop,
                     patches,
                     lambda,
                     dispersal_frac,
                     adjacency_matrix,
                     check = FALSE)

    pop <- disp$pop
    
    # migration flows: emigrants left patch i, immigrants arrived in patch j
    emigrants  <- rowSums(disp$flows) - diag(disp$flows)
    immigrants <- colSums(disp$flows) - diag(disp$flows)
    cum_emigrants  <- cum_emigrants + emigrants
    cum_immigrants <- cum_immigrants + immigrants
    
    # patches occupied
    curr_pop_size <- sapply(pop, nrow)
    occupied <- sum(curr_pop_size >= colonisation_threshold)
    
    #time for population to reach half K (carrying capapcity)
    time_halfK[is.na(time_halfK) & curr_pop_size >= half_K] <- year
    
    # arrival time and the immigrant pressure needed to reach it
    newly_colonised <- is.na(time_colonised) & curr_pop_size >= colonisation_threshold
    time_colonised[newly_colonised]        <- year
    immigrants_to_arrival[newly_colonised] <- cum_immigrants[newly_colonised]
    
    # growth rate (Eq. 5): local = reproduction only (S_t vs N_t, before dispersal);
    # net = includes migration. Undefined (NA) for empty patches
    local_growth <- ifelse(prev_pop_size > 0,
                           (post_repro_size - prev_pop_size) / prev_pop_size,
                           NA_real_)
    net_growth <- ifelse(prev_pop_size > 0,
                         (curr_pop_size - prev_pop_size) / prev_pop_size,
                         NA_real_)
    
    # genetic deficit in propagule output: realised emigrants / emigrants expected
    # with no genetic effect (m * N_t * A_t). ~1 without load, < 1 when load suppresses output
    expected_emigrants <- dispersal_frac * prev_pop_size *
      bev_holt(prev_pop_size, fecundity, carrying_capacity)
    genetic_output_ratio <- ifelse(expected_emigrants > 0,
                                   emigrants / expected_emigrants,
                                   NA_real_)
    
    # Track population statistics 
      patch_stats[[year]] <- tibble(
      year = year,
      patch = seq_along(pop),
      start_size = prev_pop_size,             # N_t
      post_repro_size = post_repro_size,      # S_t
      pop_size = curr_pop_size,               # N_{t+1} (after dispersal)
      time_half_K = time_halfK,
      g_rate = local_growth,
      g_rate_net = net_growth,
      emigrants = emigrants,
      immigrants = immigrants,
      cum_emigrants = cum_emigrants,
      cum_immigrants = cum_immigrants,
      expected_emigrants = expected_emigrants,
      genetic_output_ratio = genetic_output_ratio,
      time_colonised = time_colonised,
      immigrants_to_arrival = immigrants_to_arrival,
      patch_occupied = occupied
      )
  
  #bind population dynamics and genetics outputs 
  patch_stats_df <- bind_rows(patch_stats)
  genetic_data_df <- bind_rows(genetic_data)
  }
  
  # Return the collected data
  results <- list(
    #final_pop = pop,
    patch_stats = patch_stats_df,
    genetic_data = genetic_data_df
  )
}

