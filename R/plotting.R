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
  
  # TODO: Add plotting functions here
}

# ---------------------------------------------------------
# Plot patch size details
# ---------------------------------------------------------
plot_patch_sizes = function() {
  
  message(" > Plotting patch sizes")
  
  # Load imulation outputs: patch_sizes
  patch_sizes_dt = read_rds("compiled", "patch_sizes")
  
  # Append names of variables to column values
  plot_dt = patch_sizes_dt %>%
    mutate(across(
      .cols = c(dispersal_prob, init_frequency, lethal_effect), 
      .fns  = ~paste0(cur_column(), ": ", .x)))
  
  # Create basic plot of patch sizes over time
  g = ggplot(plot_dt) + 
    aes(x = year, 
        y = pop_size,
        colour = factor(patch), 
        linetype = factor(seed)) + 
    geom_line() + 
    facet_grid(
      cols = vars(init_frequency, lethal_effect), 
      rows = vars(dispersal_prob))
  
  # Save figure to file
  save_fig(g, "Patch sizes")
}

# ---------------------------------------------------------
# Plot allele frequency
# ---------------------------------------------------------
plot_allele_frequency = function() {
  
  message(" > Plotting allele frequency")
  
  # TODO: Plotting code here
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(g, ...) {
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    
    # Construct full figure path - including file extension
    name = paste(unlist(list(...)), collapse = " - ")
    file = paste0(o$pth$figures, name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(
      filename = file, 
      plot     = g, 
      device   = fig_format, 
      dpi      = o$save_resolution, 
      width    = o$save_width, 
      height   = o$save_height, 
      units    = o$save_units)
  }
}

