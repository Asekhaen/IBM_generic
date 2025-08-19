############################################################
# PLOTTING
#
# All plotting functionality lives here.
#
############################################################

# ---------------------------------------------------------
# Parent function for calling figure functions
# ---------------------------------------------------------
run_plotting = function() {
  
  # Only continue if specified by do_step
  if (!is.element(3, o$do_step)) return()
  
  message("* Producing figures")
  
  # Plot patch size details
  plot_patch_sizes()
  
  # Plot allele frequency
  plot_allele_frequency()
}

# ---------------------------------------------------------
# Plot patch size details
# ---------------------------------------------------------
plot_patch_sizes = function() {
  
  # ---- Load stuff ----
  
  # Load imulation outputs: patch_sizes
  patch_sizes_dt = read_rds("compiled", "patch_sizes")
  
  # ---- Create plotting dataframe ----
  
  browser()
  
  plot_dt
  
  # 
}

# ---------------------------------------------------------
# Plot allele frequency
# ---------------------------------------------------------
plot_allele_frequency = function() {
  
  # TODO: Plotting code here
}

