###########################################################
# AUXILIARY FUNCTIONS
#
# A series of helpful R functions. Note that all cluster-related
# functions now live in a separate file (cluster.R)
#
# Written by A.J.Shattock
###########################################################

# ---------------------------------------------------------
# Append columns to datatable only if they are missing
# ---------------------------------------------------------
add_if_missing = function(dt, l) {
  
  # Input should be in list format as per dplyr standard
  if (!is.list(l))
    stop("Input must be provided in list form")
  
  # Disregard any columns already existing
  l[names(l) %in% names(dt)] = NULL
  
  # Append with defined value
  for (x in names(l))
    dt[[x]] = l[[x]]
  
  return(dt)
}

# ---------------------------------------------------------
# Set as datatable and rename columns in one line
# ---------------------------------------------------------
as_named_dt = function(x, new_names) {
  
  # Convert to datatable
  dt = as.data.table(x)
  
  # Check new names are correct length
  old_names = names(dt)
  if (length(old_names) != length(new_names))
    stop("Inconsistent number of column names provided")
  
  # Set new column names
  named_dt = setnames(dt, old_names, new_names)
  
  return(named_dt)
}

# ---------------------------------------------------------
# Convert string from camel case (abcDef) to snake case (abc_def)
# ---------------------------------------------------------
camel_to_snake = function(x) {
  y = tolower(str_replace_all(x, "(?<!^)([A-Z])", "_\\1"))
  return(y)
}

# ---------------------------------------------------------
# Clear the console (cross-platform)
# ---------------------------------------------------------
clc = function() {

  # Clear a bash terminal (running R or radian)
  if (is_terminal()) system("clear")
  
  # Clear an Rstudio console
  if (!is_terminal()) cat("\014")
}

# ---------------------------------------------------------
# Clear all figures
# ---------------------------------------------------------
clf = function() {

  # Close plots
  graphics.off()

  # Close 3D rgl plots (if any)
  try(rgl::rgl.close(), silent = TRUE)
}

# ---------------------------------------------------------
# Create colour scheme
# ---------------------------------------------------------
colour_scheme = function(map, pal = NULL, n = 1, ...) {
  
  # Has colour palette been defined
  if (is.null(pal)) {
    
    # That's ok as long as it's defined within the map argument
    if (!grepl("::", map))
      stop("Palette not defined - Use 'pal = my_pal' or 'map = my_map::my_pal'")
    
    # Seperate out the map and the palette
    pal = str_remove(map, ".*\\::")
    map = str_remove(map, "\\::.*")
  }
  
  # Initiate colours variable
  colours = NULL
  
  # Built in colour schemes
  if (map == "base")
    colours = get(pal)(n, ...)	
  
  # A load of colour maps from the pals package
  #
  # See: https://www.rdocumentation.org/packages/pals/versions/1.6
  if (map == "pals")
    colours = get(pal)(n, ...)
  
  # Stylish HCL-based colour maps
  #
  # See: https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html
  if (grepl("_hcl$", map))
    colours = get(map)(palette = pal, n = n, ...)
  
  # Colour Brewer colour schemes
  if (map == "brewer")
    colours = brewer_pal(palette = first_cap(pal), ...)(n)
  
  # Viridis colour schemes
  if (map == "viridis")
    colours = viridis_pal(option = pal, ...)(n)
  
  # Throw an error if colours not yetr defined
  if (is.null(colours))
    stop("Colour map '", map, "' not recognised (supported: base, pals, hcl, brewer, viridis)")
  
  return(colours)
}

# ---------------------------------------------------------
# Reset R's most annoying default options
# ---------------------------------------------------------
default_R_options = function() {

  # All round improvements
  opts = list(
    dplyr.summarise.inform = FALSE,
    datatable.verbose = FALSE,
    stringsAsFactors = FALSE,
    scipen = 999) # Avoid scientific notation

  # Graphics settings are OS-specific
  os = get_os()

  # Windows graphics
  if (os == "windows") {
    opts$device     = "windows"
    opts$bitmapType = "windows"
  } 

  # Mac graphics
  if (os == "mac") {
    # opts$device     = "quartz"
    # opts$bitmapType = "cairo"
  }
  
  # Linux graphics
  if (os == "linux") {
    opts$device     = ifelse(capabilities()[["X11"]], "x11", "pdf")
    opts$bitmapType = "cairo"
  } 
  
  # Apply options
  do.call(options, opts)
}

# ---------------------------------------------------------
# Determine which Intercetive Development Environment (IDE) is being used
# ---------------------------------------------------------
detect_ide = function() {

  # Assume a default case
  ide = "terminal"  

  # Check for RStudio
  if (Sys.getenv("RSTUDIO") == "1")
    ide = "rstudio"
  
  # Check for VS Code
  if ("vscode" %in% Sys.getenv("TERM_PROGRAM", ""))
    ide = "vscode"

  return(ide)
}

# ---------------------------------------------------------
# Check for existance of files in a directory, ignoring subdirs
# ---------------------------------------------------------
direct_files = function(path) {

  # List immediate contents
  files = dir(
    path = path, 
    no.. = TRUE,
    full.names = TRUE, 
    recursive  = FALSE)

  # Check for any non-directories
  any_files = any(!file.info(files)$isdir)

  return(any_files)
}

# ---------------------------------------------------------
# Evaluate a string (in calling function environment) using eval
# ---------------------------------------------------------
eval_str = function(...)
  eval(parse(text = paste0(...)), envir = parent.frame(n = 1))

# ---------------------------------------------------------
# Biphasic exponential function
# ---------------------------------------------------------
exp_biphasic = function(x, peak, p, d1, d2, vmax, alpha, beta) {
  
  # Both exponential 'phases': short and long
  exp1 = exp(-x * log(2)/d1) * p
  exp2 = exp(-x * log(2)/d2) * (1-p)
  
  # Combine the phases
  bi_exp = peak * (exp1 + exp2)
  
  # Bound above and below
  bi_exp_scaled = vmax * (1 - 1 / (1 + (bi_exp / beta)^alpha))
  
  return(bi_exp_scaled)
}

# ---------------------------------------------------------
# Double exponential function
# ---------------------------------------------------------
exp_double = function(x, a, b, c, d, e) {
  y = a * exp(-b * x) + d * exp(-e * x) + c
  return(y)
}

# ---------------------------------------------------------
# Extract facets rows and columns from a ggplot object
# ---------------------------------------------------------
facet_dims = function(g) {
  g_layout = ggplot_build(g)$layout$layout
  n_rows = length(unique(g_layout$ROW))
  n_cols = length(unique(g_layout$COL))
  return(c(n_rows, n_cols))
}

# ---------------------------------------------------------
# A wrapper for tagger::tag_facets with a few extras
# ---------------------------------------------------------
facet_labels = function(g, ...) {
  
  # Default arguments
  args = list(tag = "panel", 
              tag_levels = c("A", "1"), 
              tag_suffix = "")
  
  # Overwrite these if desired
  args = list_modify(args, !!!list(...))
  
  # Call tag_facets function with these arguments
  g = g + do.call(tag_facets, args) 
  
  return(g)
}

# ---------------------------------------------------------
# More convenient logic for recode-factor process
# ---------------------------------------------------------
fct_recode = function(dt, x) {
  
  # Ensure input is in list format
  if (!is.list(x))
    stop("Input must be in list form: list(var = dict)")
  
  # Key functionality - recode then set as factors
  x1 = "mutate(X = recode(X, !!!dict))"
  x2 = "mutate(X = factor(X, dict))"
  
  # Iterate through variables
  vars = names(x)
  for (var in vars) {
    
    # Extract dictionary for this variable
    dict = x[[var]]
    
    # Replace placeholder with variable name
    s1 = str_replace_all(x1, "X", var)
    s2 = str_replace_all(x2, "X", var)
    
    # Evaluate string to update dt
    eval_str("dt %<>%", s1, "%>%", s2)
  }
  
  return(dt)
}

# ---------------------------------------------------------
# Platform specific file separator - for readability
# ---------------------------------------------------------
file_sep = function() {
  platform_file_sep = .Platform$file.sep
  return(platform_file_sep)
}

# ---------------------------------------------------------
# Size of all files in a directory, returns a vector
# ---------------------------------------------------------
file_sizes = function(pth) {
  
  # NOTE: Only works on UNIX systems
  
  # Construct linux find command
  command = paste("find", shQuote(pth), "-type f -printf '%s\n'")
  
  # Exectute shell command
  file_sizes = as.numeric(system(command, intern = TRUE))
  
  return(file_sizes)
}

# ---------------------------------------------------------
# Capitalise first letter of a string (or vector of strings)
# ---------------------------------------------------------
first_cap = function(string) {
  string_cap = paste0(toupper(substring(string, 1, 1)), substring(string, 2))
  return(string_cap)
}

# ---------------------------------------------------------
# Format heterogeneous styles of dates
# ---------------------------------------------------------
format_date = function(dates, convert = "ymd") {
  styles = c("dmy", "dmY", "ymd", "Ymd")
  dates = parse_date_time(dates, styles)
  dates = get(convert)(dates)
  return(dates)
}

# ---------------------------------------------------------
# Identify operating system
# ---------------------------------------------------------
get_os =  function() {

  # Define dictionary mapping
  dict_os = c(
    "Windows" = "windows",
    "Linux"   = "linux",
    "Darwin"  = "mac")

  # Load name from system info
  sys_name = Sys.info()[["sysname"]]

  # Match to dictionary
  os_name = ifelse(
    test = sys_name %in% names(dict_os), 
    yes  = dict_os[sys_name], 
    no   = "unknown")

  return(os_name)
}

# ---------------------------------------------------------
# Interpolate time series trends
# ---------------------------------------------------------
interp_ts_trend = function(dt) {
  interp_dt = dt %>%
    model(lm = TSLM(log(value) ~ trend())) %>%
    interpolate(dt)
  return(interp_dt)
}

# ---------------------------------------------------------
# Check whether running in Rstudio or bash terminal
# ---------------------------------------------------------
is_terminal = function() {
  
  # Check for RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (rstudioapi::isAvailable()) {
      
      # Return false
      return(FALSE)
    }
  }
  
  # Check for TTY properties
  tty = suppressWarnings(system("tty", intern = TRUE))
  
  # Check properties of terminal
  if (is.character(tty)) {
    if (file.exists(tty)) {
      
      # Return true
      return(TRUE)
    }
  }
  
  # Otherwise unknown
  return(NULL)
}

# ---------------------------------------------------------
# Convert list to datatable
# ---------------------------------------------------------
list2dt = function(x, ...) {
  dt = rbindlist(lapply(x, as.data.table), ...)
  return(dt)
}

# ---------------------------------------------------------
# General log-likelihood for a Gaussian distribution
# ---------------------------------------------------------
ll_gaussian = function(data, model, sigma, weight) {
  
  # Calculate residuals
  residuals = data - model
  
  # Apply time weights to residuals (for periods of strong seasonality)
  weighted_residuals = residuals * weight
  
  # Calculate log-likelihood assuming Gaussian distributed errors
  ll = dnorm(
    x    = weighted_residuals, 
    mean = 0, 
    sd   = sigma, 
    log  = TRUE)
  
  # Sum the log likelihood
  ll_total = -sum(ll)
  
  return(ll_total)
}

# ---------------------------------------------------------
# Logistic function
# ---------------------------------------------------------
logistic = function(x, slope, mid, lower = 0, upper = 1) {
  y = upper + (lower - upper) / (1 + (x / mid) ^ slope)
  return(y)
}

# ---------------------------------------------------------
# Inverse logistic function (solving for x)
# ---------------------------------------------------------
logistic_inv = function(y, slope, mid, lower, upper) {
  x = mid * ((lower - y) / (y - upper)) ^ (1 / slope)
  return(x)
}

# ---------------------------------------------------------
# Reverse logistic function
# ---------------------------------------------------------
logistic_rev = function(x, slope, mid, lb, ub) {
  y = 1 - logistic(x, slope, mid, lower = 1 - ub, upper = 1 - lb)
  return(y)
}

# ---------------------------------------------------------
# Parameter transformation: put a probability on the real number line
# ---------------------------------------------------------
logit = function(p) {
  z = log(p / (1 - p))
  return(z)
}

# ---------------------------------------------------------
# Parameter transformation: inverse of the above
# ---------------------------------------------------------
logit_inv = function(p_logit) {
  z = exp(p_logit) / (exp(p_logit) + 1)
  return(z)
}

# ---------------------------------------------------------
# Moving average that doesn't truncate
# ---------------------------------------------------------
moving_avg = function(x, source) {
  
  # Uses caTools package
  y = runmean(
    x = x,
    k = o$data_ma[[source]],
    endrule = "keep")
  
  return(y)
}

# ---------------------------------------------------------
# Simple wrapper for number of unique observations
# ---------------------------------------------------------
n_unique = function(x) length(unique(x))

# ---------------------------------------------------------
# Wrapper for lapply that also extracts element name
# ---------------------------------------------------------
napply = function(x, fn, ...) {
  y = lapply(seq_along(x), function(i) fn(x[[i]], name = names(x)[i], ...))
  return(y)
}

# ---------------------------------------------------------
# Normalise a vector of values to between 0 and 1
# ---------------------------------------------------------
normalise_0to1 = function(x, x_min = NULL, x_max = NULL, direction = "forward") {
  
  if (!tolower(direction) %in% c("forward", "backward"))
    stop("Normalisation direction must be either 'forward' or 'backward'")
  
  # Forward normalisation
  if (tolower(direction) == "forward") {
    
    # Take bounds from data unless given
    if (is.null(x_min)) x_min = min(x)
    if (is.null(x_max)) x_max = max(x)
    
    # Normalisation equation
    y = (x - x_min) / (x_max - x_min)
    
    # Append original min and max values, needed to backtransform
    attributes(y)$x_min = x_min
    attributes(y)$x_max = x_max
  }
  
  # Backward normalisation
  if (tolower(direction) == "backward") {
    
    # Take bounds from attitubutes of pre-normalised data unless given
    if (is.null(x_min)) x_min = attributes(x)$x_min
    if (is.null(x_max)) x_max = attributes(x)$x_max
    
    # Rearrange equation to solve for x
    #
    # NOTE: as.vector removes all attributes
    y = as.vector(x * (x_max - x_min) + x_min)
  }
  
  return(y)
}

# ---------------------------------------------------------
# Equivalent of paste, but with an underscore instead of space
# ---------------------------------------------------------
paste1 = function(...) paste(..., sep = "_")

# ---------------------------------------------------------
# Sub a directory name within a file path string
# ---------------------------------------------------------
pth_replace = function(pth, dir, dir_replace, sep = file_sep()) {
  new_pth = gsub(paste0(sep, dir, sep), 
                 paste0(sep, dir_replace, sep), pth)
  return(new_pth)
}

# ---------------------------------------------------------
# Suppress output from a function call
# ---------------------------------------------------------
quiet = function(x) { 
  sink_con = file("sink.txt")
  sink(sink_con, type = "output")
  sink(sink_con, type = "message")
  on.exit(sink(type   = "output"))
  on.exit(sink(type   = "message"), add = TRUE)
  on.exit(file.remove("sink.txt"),  add = TRUE)
  invisible(force(x)) 
}

# ---------------------------------------------------------
# Convenience wrapper for readRDS
# ---------------------------------------------------------
read_rds = function(pth, ..., err = TRUE) {
  
  # Special use case: pth is the full .rds file path
  if (grepl(".*\\.rds$", pth)) {
    full_path = pth
    
  } else { # Otherwise standard use case
    
    # Construct path and file name using inputs
    file_path = o$pth[[pth]]
    file_name = paste(unlist(list(...)), collapse = "_")
    
    # Concatenate full .rds file path
    full_path = paste0(file_path, file_name, ".rds")
  }
  
  # Check whether file exists
  exists = file.exists(full_path)
  
  # If file exists, load it
  if (exists)
    x = readRDS(file = full_path)
  
  # If file does not exist
  if (!exists) {
    
    # Construct error / warning message
    err_message = paste("Unable to load file '", full_path, "'")
    
    # Either throw an error or warning depending on err argument
    if (err) stop(err_message)
    if (!err) warning(err_message)
    
    # Return out trivial result
    x = NULL
  }
  
  return(x)
}

# ---------------------------------------------------------
# Simple wrapper for read_yaml -> list2dt
# ---------------------------------------------------------
read_yaml_dt = function(x) list2dt(read_yaml(x))

# ---------------------------------------------------------
# Simple wrapper for read_yaml_dt from input dir
# ---------------------------------------------------------
read_yaml_input = function(x, dt = TRUE) {
  
  # Append file name with input path
  path = paste0(o$pth$input, x, ".yaml")
  
  # Call regular read_yaml function
  if (dt)  y = read_yaml_dt(path)
  if (!dt) y = read_yaml(path)
  
  return(y)
}

# ---------------------------------------------------------
# Wrapper for consistent behaviour of base::sample when length(x) is one
# ---------------------------------------------------------
sample_vec = function(x, ...) x[sample.int(length(x), ...)]

# ---------------------------------------------------------
# Convenience wrapper for saveRDS
# ---------------------------------------------------------
save_rds = function(x, pth, ...) {
  
  # Special use case: pth is the full .rds file path
  if (grepl(".*\\.rds$", pth)) {
    full_path = pth
    
  } else { # Otherwise standard use case
    
    # Construct path and file name using inputs
    file_path = o$pth[[pth]]
    file_name = paste(unlist(list(...)), collapse = "_")
    
    # Concatenate full .rds file path
    full_path = paste0(file_path, file_name, ".rds")
  }
  
  # Save as an RDS
  saveRDS(x, file = full_path)
}

# ---------------------------------------------------------
# Initiate progress bar with normal-use options
# ---------------------------------------------------------
start_progress_bar = function(n = 100) {

  # Min and max size of progress bar
  min_size = 40
  max_size = 70

  # Get current terminal size
  term_size = as.numeric(system("tput cols", intern = TRUE))
  
  # Set progress bar size, forcing between min and max
  size = min(max(term_size, min_size), max_size)
  
  # Initiate progress bar from progress package
  pb = progress_bar$new(
    format     = " [:bar] :percent (:eta)",
    total      = n,     # Number of tasks to complete
    complete   = "-",   # Completion bar character
    incomplete = " ",   # Incomplete bar character
    current    = ">",   # Current bar character
    clear      = TRUE,  # If TRUE, clears the bar when finished
    width      = size)   # Width of the progress bar
  
  return(pb)
}

# ---------------------------------------------------------
# Convert comma-seperated string to vector of elements
# ---------------------------------------------------------
str2vec = function(x, v) {
  x[[v]] = x[[v]] %>% 
    str_split(",", simplify = TRUE) %>% 
    str_remove_all(" ")
  return(x)
}

# ---------------------------------------------------------
# Bi-directional setdiff - elements not in both x and y
# ---------------------------------------------------------
symdiff = function(x, y) setdiff(union(x, y), intersect(x, y))

# ---------------------------------------------------------
# Format a number with thousand mark separators
# ---------------------------------------------------------
thou_sep = function(val) {
  format_val = format(val, scientific = FALSE,
                      trim = TRUE, 
                      drop0trailing = TRUE, 
                      big.mark = ",")
  return(format_val)
}

# ---------------------------------------------------------
# Format time taken into h-m-s string
# ---------------------------------------------------------
time_taken = function(id) {
  
  # Time since last tic
  time = toc(quiet = TRUE)
  
  # Convert to lubridate duration
  duration = seconds_to_period(time$toc - time$tic)
  
  # Construct a human-readable string
  string = sprintf(
    fmt = "%02dh %02dm %02ds",
    duration@hour, 
    duration@minute, 
    floor(duration@.Data))
  
  return(string)
}

# ---------------------------------------------------------
# Load an file if it exists, throw an error if not
# ---------------------------------------------------------
try_load = function(pth, file, msg = NULL, type = "rds", throw_error = TRUE, sep = FALSE) {
  
  # Initiate trivial output
  file_contents = NULL
  
  # Set default error message
  if (is.null(msg))
    msg = "Cannot load file"
  
  # Switch case for loading function
  loading_fnc = switch(
    tolower(type), 
    
    # Support both RDS and CSV
    "rds" = "readRDS", 
    "csv" = "read.csv",
    
    # Throw an error if anything else requested
    stop("File type '", type, "' not supported")
  )
  
  # Concatenate path and file name
  file_name = paste0(pth, ifelse(sep, file_sep(), ""), file, ".", type)
  
  # If file doesn't exist, throw an error if desired
  if (!file.exists(file_name) && throw_error == TRUE)
    stop(msg, " [missing: ", file_name, "]")
  
  # If file exists, try to load it
  if (file.exists(file_name)) {
    
    # Get the loading function and attempt to load file
    file_contents = tryCatch(
      get(loading_fnc)(file_name),
      
      # Catch the error - we may not want to throw it
      error = function(e) {
        
        # Throw descriptive error if desired
        if (throw_error == TRUE) 
          stop(msg, " [unreadable: ", file_name, "]")
      }
    )
  }
  
  return(file_contents)
}
  
# ---------------------------------------------------------
# Adapted Weibull function - as used in OpenMalaria
# ---------------------------------------------------------
weibull_adapted = function(t, a, b, c) {
  y = a * exp(-(t/b)^c * log(2))
  return(y)
}

