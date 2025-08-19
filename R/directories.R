############################################################
# SET DIRECTORIES
#
# Set and get directories in one place in the name of consistency
# and ease. Also creates any directories that do not currently exist.
#
# Outputs a list of relevant directories (within o$pth) which 
# can be referenced elsewhere.
#
############################################################

# ---------------------------------------------------------
# Define paths for project inputs and outputs
# ---------------------------------------------------------
set_dirs = function(o) {
  
  # Initiate file path lists
  pth = list()
  
  # ---- Code and resource locations ----
  
  # Current directory - code library
  pth$code = getwd()
  
  # Directory for input files, snippets, and cluster logs
  pth$input = file.path(dirname(pth$code), "input")
  pth$log   = file.path(dirname(pth$code), "log")
  
  # Path for all output files (linked to output version)
  pth$output = set_output_dir(o, pth)
  
  # ---- Output directories ----
  
  # Paths to compiled data for use in history matching and fitting
  pth$compiled = file.path(pth$output, "0_compiled")
  
  # Paths to post-processed files
  pth$sims = file.path(pth$output, "1_simulations*")
  
  # Path to analysis files
  pth$analysis = file.path(pth$output, "2_analysis")
  
  # Path to results and figures
  pth$figures = file.path(pth$output, "3_figures")
  
  # ---- Create directories if needed ----

  # Setnix only: Set up for directories to be stored on acacia
  pth = identify_acacia_dirs(pth)  # See acacia.R
  
  # Make all output directories
  pth = make_output_dirs(pth)
  
  # Append paths to o list
  o$pth = pth
  
  return(o)
}

# ---------------------------------------------------------
# Set path to output directory
# ---------------------------------------------------------
set_output_dir = function(o, pth) {
  
  # Output version (linked to output directory name)
  dir_name = paste1("output", o$output_version)
  
  # Path to output pth: local or scicore - same as code repo (keeps things tidy)
  if (o$hpc$name != "setonix")
    dir_path = dirname(pth$code)
  
  # Path to output pth: setonix
  if (o$hpc$name == "setonix") {

    # NOTE: Different set up for Pawsey - code repo in software, output in scratch

    # Name of this repo
    repo = basename(dirname(pth$code))

    # Other details for constructuing path
    proj = o$hpc$project
    user = o$user

    # Construct path to scratch where output should be stored
    dir_path = file.path("", "scratch", proj, user, repo)
  }
  
  # Append output pth name to for full path
  pth_output = file.path(dir_path, dir_name)
  
  return(pth_output)
}

# ---------------------------------------------------------
# Make all output directories if they do not already exist
# ---------------------------------------------------------
make_output_dirs = function(pth) {
  
  # Reset variable
  x = pth
  
  # Remove storage and ignore directories for this process
  x = x[!names(x) %in% c("storage", "ignore")]
  
  # Iterate through dirs
  for (i in names(x)) {
    p = x[[i]]
    
    # If it does not already exist, create it
    if (!dir.exists(p))
      dir.create(p, recursive = TRUE)
    
    # Add a file separator to end of dir path
    x[[i]] = paste0(p, file_sep())
  }

  # Reappend strorage and ignore fields
  x$storage = pth$storage
  x$ignore  = pth$ignore
  
  return(x)
}

