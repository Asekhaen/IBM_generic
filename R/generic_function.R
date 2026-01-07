


# =================================
# POPULATION INITIALISATION 
# =================================

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


# ===============================================
# REPRODUCTION AND GENETIC INHERITANCE  
# ===============================================

# captures density-dependent reproduction using Beverton-Holt model, and genetic 
# (with linkage) inheritance 

growth <- function(pop_patches, 
                   n_loci,
                   carrying_capacity,
                   fecundity,
                   lethal_effect,
                   complete_sterile,
                   linkage,
                   sim_years) {
  #browser()
  #if(sim_years == 4) browser()
  updated_pop_patches <- list()
  #updated_fitness <- list ()
  for (i in seq_along(pop_patches)) {
    pop <- pop_patches[[i]]  
    
    # reproduction  
    n.pop <- nrow(pop)
    
    if (n.pop > 0){
      
      
      # exp_fecundity <- fec_dd(n.pop, dd_rate, prob_survival)
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
        
        homozygous <- rowSums(
          ((pop$allele1 + pop$allele2) == 2) |
            ((pop$mate_allele1 + pop$mate_allele2) == 2)        # if either parent is homozygous in at least one locus
        ) > 0
        
        sterile <- as.numeric(!homozygous)
        n_offspring <- act_fecundity * sterile
        

      } else {
        n_offspring <- act_fecundity
      }
      
    }   else {
      # If not, set offspring count to 0
      n_offspring <- rep(0, n.pop)
    }
    
    
     total_offspring <- sum(n_offspring)
     
    

 # GENETIC INHERITANCE

     
# function with option to choose between inheritance type (linkage and without linkage)

     if (total_offspring > 0){

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

         # Offspring tibble
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

     # if turned on, this condition statement simulates lethal effect for individuals with 
     # homologous deleterious alleles (in at least one of loci of either parents)
     
     if (lethal_effect){
       homozygous_lethal <- (pop$allele1 == 1) & (pop$allele2 == 1)
       any_homozygous <- rowSums(homozygous_lethal) > 0
       pop <- pop[!any_homozygous,]
     }

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
  
  dispersed_pop
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
      
      
      tibble::tibble(
        patch = patch_id,
        year  = year,
        locus = 1:n_loci,
        AA = AA_count,
        Aa = Aa_count,
        aa = aa_count,
        freq_A = freq_A,
        freq_a = freq_a
      )
    })
    

    #Dispersal
    pop <- dispersal(pop,
                     patches,
                     lambda,
                     dispersal_frac,
                     adjacency_matrix,
                     check = FALSE)

    # patches occupied
    curr_pop_size <- sapply(pop, nrow)
    occupied <- sum(curr_pop_size >= establish_threshold)
    
    #time for population to reach half K (carrying capapcity)
    time_halfK[is.na(time_halfK) & curr_pop_size >= half_K] <- year
    
    
    # growth rate 
    growth_rate <- ifelse(prev_pop_size > 0,
                          ((curr_pop_size - prev_pop_size) / prev_pop_size),
                          0)
    
    # Track population statistics (using population statistics before dispersal? how about post dispersal?)
    patch_stats[[year]] <- tibble(
      year = year,
      patch = seq_along(pop),
      pop_size = curr_pop_size,
      time_half_K = time_halfK,
      g_rate =  growth_rate,
      patch_occupied = occupied
    )
  

  #bind population dynamics outputs 
  patch_stats_df <- bind_rows(patch_stats)
  # estimate spread rate or speed
  patch_stats_df <- invasion_speed(patch_stats_df)
  

  #genotype and allele frequency
  genetic_data_df <- bind_rows(genetic_data)
  }
  
  
  #-------------------------------------
  # Plot landscape (igraph package)
  #-------------------------------------

# if (adjacency_matrix) {
#   # Plot One-dimensional landscape with nearest-neighbour
# 
#   x <- 1:patches
#   y <- rep(0, patches)
# 
# 
#   #pdf("1D_landscape.pdf", width = 8, height = 6)
#   png("1D_landscape.png", width = 1200, height = 800, res = 150)
# 
#   plot(
#     x, y,
#     pch = 21, bg = "skyblue",
#     cex = 6,
#     xlab = "Patches",
#     ylab = "",
#     yaxt = "n",
#     xlim = c(0, patches+1),
#     ylim = c(-0.5, 0.5),
#     main = "1D stepping-stone"
#   )
# 
#   # Label patches
#   text(x, y, labels = paste0("P", 1:patches), cex = 1, col = "black")
# 
#   dev.off()
# } else {
#   set.seed(26)
#   #plot metapopulation network landscape
# 
#   coords <- create_coordinates(patches)
#   # dispersal matrix
#   dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
#   #exponential dispersal kernel
#   dispersal_kernel <- exp(-lambda * dist_matrix)
#   # set the diagonal elements to 0 to prevent self-dispersal
#   diag(dispersal_kernel) <- 0
# 
# 
#   #Plot using weights for edge width (using R package "igraph")
#   g <- graph_from_adjacency_matrix(dispersal_kernel, mode="undirected", weighted=TRUE, diag=FALSE)
# 
#   png("metapop_landscape.png", width = 1200, height = 800, res = 150)
# 
#   plot(
#     g,
#     layout = coords,
#     vertex.size = 35,
#     vertex.color = "skyblue",
#     vertex.label = paste0("P", 1:patches),
#     vertex.label.color = "black",
#     edge.width = E(g)$weight * 10,  # scale for visibility
#     edge.color = rgb(0,0,0, alpha=E(g)$weight),
#     main = "Metapopulation network"
#   )
# dev.off()
# }

  # Return the collected data
  results <- list(
    #final_pop = pop,
    patch_stats = patch_stats_df,
    genetic_data = genetic_data_df
  )
}

