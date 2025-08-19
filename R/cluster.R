###########################################################
# CLUSTER FUNCTIONS
#
# All cluster-related R functions in one place.
#
# Written by A.J.Shattock
###########################################################

# ---------------------------------------------------------
# Pull a few details about the user and their setup
# ---------------------------------------------------------
identify_user = function() {
  
  # Initiate list to store these details
  o = list()
  
  # User-specific details (for checking cluster jobs)
  o$user = Sys.info()[["user"]]
  
  # Identify High Performance Computing system
  o$hpc = identify_hpc()
  
  # Number of cores available (for parallel loading)
  o$cores = detectCores()
  
  return(o)
}

# ---------------------------------------------------------
# Identify which High Performance Computing system we're using
# ---------------------------------------------------------
identify_hpc = function() {
  
  # We do this by assessing module paths
  str = system("echo $MODULEPATH", intern = TRUE)
  
  # Search for 'scicore' or 'setonix'
  hpc = list(
    name = dplyr::case_when(
      grepl("setonix", str) ~ "setonix",
      TRUE ~ "none"))
  
  # Also append node name (to identify login or compute node)
  if (hpc$name != "none")
    hpc$node = system("echo $(hostname)", intern = TRUE)
  
  # If running on Setonix, extract project ID
  if (hpc$name == "setonix")
    hpc$project = sub(".*/(pawsey[0-9]+)/.*", "\\1", str)
  
  return(hpc)
}

# ---------------------------------------------------------
# Submit jobs to the cluster and wait until they are complete
# ---------------------------------------------------------
submit_cluster_jobs = function(n_jobs, bash_file, ...) {
  
  # Skip if no jobs to run
  if (n_jobs > 0) {
    
    # Throw an error if no HPC identified
    if (o$hpc$name == "none")
      stop("You do not appear to be connected to a HPC - aborting cluster submission")
    
    # Extract name of function to call
    job_name = list(...)[[1]]
    
    # ---- Prepare for cluster submission ----
    
    # Check if user is currently running any cluster jobs (see auxiliary.R)
    check_cluster_jobs(action = o$cluster_conflict_action)
    
    # Create a new log file for the cluster jobs
    log_file = create_bash_log(
      pth = o$pth$log, 
      log = o$log_file, 
      err = o$err_file)
    
    # ---- Construct and execute array job submission ----
    
    # Construct SBATCH array and options strings
    sbatch_array   = set_sbatch_array(n_jobs)
    sbatch_options = set_sbatch_options(job_name)
    
    # Construct list of all inputs to bash file
    bash_inputs = set_bash_inputs(n_jobs, job_name, log_file)
    
    # Concatenate system command
    sys_command = paste("sbatch", sbatch_options, sbatch_array, bash_file, bash_inputs)
    
    # message(str_replace_all(sys_command, " ", "\n"))
    
    # Invoke this command
    sys_status = system(sys_command, intern = TRUE)
    
    # ---- Create report upon completion ----
    
    # Extract array ID
    array_id = str_extract(sys_status, "[0-9]+")
    
    message(" > Submitting jobs (array ID: ", array_id, ")")
    
    # Path to report file written when all job complete
    report_file = paste0(o$pth$log, "report_", array_id, ".txt")
    
    # Basic HPC-specific settings for running job (no job size details needed)
    report_options = set_sbatch_options(job_name, job_size = FALSE)
    
    # Construct command for creating report upon completion of array job
    report_after   = paste0("--dependency=afterany:", array_id)
    report_command = paste("sbatch", report_after, report_options,
                           "report.sh", array_id, report_file)
    
    # Invoke this command
    report_status = system(report_command, intern = TRUE)
    
    # Wait for all cluster jobs to complete
    wait_for_report(report_file, n_jobs, log_file)
    
    # Print cluster report to console
    cluster_report(report_file)
  }
}

# ---------------------------------------------------------
# Check if user is currently running any cluster jobs
# ---------------------------------------------------------
check_cluster_jobs = function(action = "error") {
  
  # Check number of running and pending jobs
  n_jobs = n_cluster_jobs(user = o$user)
  
  # If this is non-zero then action needed
  if (sum(unlist(n_jobs)) > 0) {
    
    # Throw an error
    if (action != "none")
      stop("You currently have jobs submitted to the cluster, this may lead to unexpected results.")
  }
}

# ---------------------------------------------------------
# Number of running and pending jobs on the cluster
# ---------------------------------------------------------
n_cluster_jobs = function(user) {
  
  # Base sq command for user
  sq = paste("squeue -u", user)
  
  # Concatenate commands for running and pending jobs
  slurm_running  = paste(sq, "-t running | wc -l")
  slurm_pending  = paste(sq, "-t pending | wc -l")
  
  # Interactive jobs depend on HPC system
  slurm_ondemand = list(
    setonix = paste(sq, "-n interactive | wc -l"),
    scicore = paste(sq, "-q interactive | wc -l"))
  
  # Function to get number of jobs (minus 1 to remove header row)
  get_jobs_fn = function(x) as.numeric(system(x, intern = TRUE)) - 1
  
  # System call to determine number of slurm processes
  n_running  = get_jobs_fn(slurm_running)
  n_pending  = get_jobs_fn(slurm_pending)
  n_ondemand = get_jobs_fn(slurm_ondemand[[o$hpc$name]])  # Interactive jobs
  
  # Compile into list
  #
  # NOTE: ondemand jobs are considered 'running', so discount these
  n_jobs = list(running = n_running - n_ondemand, pending = n_pending)
  
  return(n_jobs)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs)
# ---------------------------------------------------------
create_bash_log = function(pth, log = NULL, err = NULL) {
  
  # Repeat process for log and error files
  for (this_name in c(log, err)) {
    this_file = paste0(pth, this_name)
    
    # Delete any existing file
    if (file.exists(this_file))
      file.remove(this_file)
    
    # Pause for file system to catch up 
    Sys.sleep(0.5)
    
    # Create a new file with a single header row
    cat("id\n", file = this_file)
    
    # Pause for file system to catch up 
    Sys.sleep(0.5)
  }
  
  # Concatenate log path and file name
  log_file = paste0(pth, log)
  
  return(log_file)
}

# ---------------------------------------------------------
# Construct sbatch array command for running in parallel
# ---------------------------------------------------------
set_sbatch_array = function(n_jobs) {
  
  # Shorthand for number of jobs that can be done on a single node
  n_tasks = o$job_size$ntasks
  
  # Number of nodes actually needed
  n_nodes = min(ceiling(n_jobs / n_tasks), o$job_limit)
  
  # This defines the size of the array we request
  sbatch_array = paste0("--array=1-", n_nodes)
  
  message("  - Array jobs: ", thou_sep(n_nodes),
          " (of ", thou_sep(n_tasks), " tasks each)")
  
  return(sbatch_array)
}

# ---------------------------------------------------------
# Format user-defined options into sbatch-interpretable options
# ---------------------------------------------------------
set_sbatch_options = function(job_name, job_size = TRUE) {
  
  # Get sbatch options for running on this HPC
  sbatch = o$job_settings[[o$hpc$name]]
  
  # Append job size options by default
  if (job_size == TRUE)
    sbatch = c(sbatch, o$job_size)
  
  # Some options are specific to the type of job (eg time and memory)
  job_specific = names(sbatch[map(sbatch, ~length(.x)) > 1])
  
  # Select the appropriate value for this type of job
  for (i in job_specific)
    sbatch[i] = sbatch[[i]][[job_name]]
  
  # Convert to a single string of multiple sbatch options
  sbatch_options = sbatch %>%
    map(~as.character(.x)) %>%
    enframe() %>%
    unnest(value) %>%
    mutate(name = paste0("--", name)) %>%
    unite("x", sep = "=") %>%
    pull() %>%
    paste(collapse = " ")
  
  return(sbatch_options)
}

# ---------------------------------------------------------
# Construct list of all inputs to bash file
# ---------------------------------------------------------
set_bash_inputs = function(n_jobs, job_name, log_file) {
  
  # Easy access key cluster variables
  hpc_name = o$hpc$name
  n_cores  = o$job_size$ntasks
  
  # R version being run
  r_version = paste0(R.Version()$major, ".", R.Version()$minor)
  
  # Extract and format arguments to bash file (including R version & log file)
  bash_inputs = paste(
    job_name,   # Input 1: Job type (eg history, fitting, future)
    n_jobs,     # Input 2: Total number of jobs to run
    n_cores,    # Input 3: Number of cores to operate per node
    hpc_name,   # Input 4: High Performance Computing system we're using
    r_version,  # Input 5: R version to load and run
    log_file)   # Input 6: Log file to write to when job complete
  
  return(bash_inputs)
}

# ---------------------------------------------------------
# Wait until all cluster jobs have finished
# ---------------------------------------------------------
wait_for_report = function(report_file, n_jobs, log_file, wait = 1) {
  
  # Wait for log file to be created
  while (!file.exists(log_file)) 
    Sys.sleep(wait)
  
  # Initiate a progress bar
  pb = start_progress_bar(100)
  pb$update(0)
  
  # Initiate flag that report file exists, and job counter
  x = FALSE
  k = 0
  
  # Wait for all jobs to write to log file
  while (x == FALSE) {
    
    # Check whether report file exists
    x = file.exists(report_file)
    
    # Read number of lines in log file
    k = nrow(read.table(log_file)) - 1
    
    # NOTE: k is is the number of completed jobs (sucessfully
    #       or not), but doesn't capture jobs that may have
    #       timed out or run out of memory (as per slurm options)
    
    # Update progress bar safely (only if it's not finished)
    if (!pb$finished) 
      pb$update(k / n_jobs)
    
    # Wait before testing again
    Sys.sleep(wait)
  }
  
  # Terminate progress bar if necessary
  if (!pb$finished) 
    pb$terminate()
}

# ---------------------------------------------------------
# Display details of clutser report
# ---------------------------------------------------------
cluster_report = function(report_file) {
  
  message(" > Preparing job report")
  
  # TODO: Rewrite this. If task info not available, look instead at cluster_log.txt
  
  # Load report from file (sacct form)
  report_raw = read.table(report_file, header = TRUE) 
  
  # Analyse the file, and count each exit status
  report_vct = report_raw %>%
    rename_with(tolower) %>%
    filter(!grepl("---", jobid)) %>%
    count(state) %>%
    arrange(-n) %>%
    mutate(p = paste0(round(100 * n / sum(n), 2), "%)"), 
           n = paste0(": ", thou_sep(n), " ("), 
           s = paste0(first_cap(tolower(state)), n, p)) %>%
    pull(s)
  
  # Combine into a string to be reported back to user
  report_str = paste0("  - ", report_vct) %>%
    paste(collapse = "\n")
  
  # Report back to user
  message(report_str)
}

# ---------------------------------------------------------
# Check for cluster errors, stop if any found
# ---------------------------------------------------------
stop_if_errors = function(pth, err, err_tol = 0, msg = NULL) {
  
  # Set default error message
  if (is.null(msg))
    msg = "Fatal errors when running cluster jobs"
  
  # Check the error file exists
  err_file = paste0(pth, err)
  if (file.exists(err_file)) {
    
    # Load errors from file and convert to vector
    errors = readLines(err_file)[-1]
    
    # Stop if any errors, and display them all
    if (length(errors) > err_tol)
      stop(msg, " (", length(errors), " errors):\n", 
           paste(errors, collapse = "\n"))
  }
}

