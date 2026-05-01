

# 1000 loci
n_loci = 1000
n_load <- c(0.001, 0.0025, 0.005, 0.01, 0.25)
# n_load <- c(0.05, 0.1, 0.25, 0.5,  0.75)
n_freq <- c(0.001000250, 0.001582127, 0.002238868, 0.003170219, 0.016959973)




# 100 loci
n_loci = 100
n_load <- c(0.001, 0.0025, 0.005, 0.01, 0.25)
n_freq <- c(0.003163061, 0.005003098, 0.007079842, 0.010024884, 0.053597450)



#function to cal del allele frequency from load
calc_q <- function(n_load, n_loci) {
  sqrt(1 - (1 - n_load)^(1 / n_loci))
}

# calc_q(0.4, 200)


#function to cal load from del allele frequency
genetic_load <- function(n_freq, n_loci){
    L <- 1 - ((1 - n_freq^2)^n_loci)
  return(L)
}

# genetic_load(n_freq, n_loci)
