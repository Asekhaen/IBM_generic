###########################################################
# ACACIA FUNCTIONS
#
# All functions related to Acacia storage in one place. 
# Write, read, and diagnostic functionality. Only relevant
# if running on Setonix.
#
# Written by A.J.Shattock
###########################################################

# ---------------------------------------------------------
# Parent function to upload/download content to/from Acacia storage
# ---------------------------------------------------------
run_acacia = function(direction) {
  
  # Only continue if specified by do_step
  if (!is.element(5, o$do_step)) return()
  
  # Ensure valid direction
  if (!direction %in% c("upload", "download"))
    stop("Input 'direction' must be either 'upload' or 'download'")
  
  # Either upload of download results from Acacia
  get(paste1("acacia", direction))()
}

# ---------------------------------------------------------
# Primary function: upload files to Acacia storage
# ---------------------------------------------------------
acacia_upload = function() {
  
  # If not running on Setonix, return out
  if (o$hpc$name != "setonix") {
    
    message("! Must be connected to Setonix to upload to Acacia !")
    
    return()
  }
  
  message("* Uploading files to Acacia")
  
  # Create a bucket on Acacia (if it doesn't already exist)
  create_acacia_link()
  
  # ---- Determine files to copy ----
  
  # All existing output directories (for this age group)
  dirs_all = dir_ls(
    path = o$pth$output, 
    type = "directory", 
    recurse = TRUE)
  
  # Paths to ignore when copying
  dirs_ignore = o$pth$ignore %>%
    map(~o$pth[[.x]]) %>%
    unlist() %>%
    str_remove("/$")
  
  # Directories to be copied (excl those with no direct files)
  dirs_copy = setdiff(dirs_all, dirs_ignore) %>% 
    keep(direct_files)
  
  # List all files in each of these directories
  files_copy = dirs_copy %>%
    map(~dir_ls(.x, type = "file", recurse = FALSE)) %>%
    unlist() %>%
    str_remove(o$pth$output)
  
  # ---- Intermediate file ----
  
  # Path to intermediate file
  path_copy = paste0(o$pth$log, "acacia.txt")
  
  # Intermediate file contains list of files to uploaded
  writeLines(files_copy, path_copy)
  
  # TODO: Throw warning (or error) if trying to write too many files
  
  # ---- Copy files ----
  
  # Construct rclone command to upload unchanged files
  rclone_cmd = paste(
    "rclone copy",
    o$pth$output, 
    get_bucket_path(),  
    "--files-from", path_copy, 
    "--progress", 
    "--checksum")
  
  # Execute rclone command
  system(rclone_cmd, wait = TRUE)
  
  # Remove intermediate file
  invisible(file.remove(path_copy))
}

# ---------------------------------------------------------
# Primary function: download files from Acacia storage
# ---------------------------------------------------------
acacia_download = function() {
  
  # NOTE: This function isn't robust to any difference in username
  #       on local vs HPC. In my example (ashattock), they are the
  #       same. If different, a data dict will be needed to convert.
  
  message("* Downloading files from Acacia")
  
  # Throw an error if rclone is not installed
  if (Sys.which("rclone") == "")
    stop("Must have rclone installed locally")
  
  # Build the rclone command
  rclone_cmd = paste(
    "rclone copy", 
    get_bucket_path(),
    o$pth$output, 
    "--local-no-set-modtime", 
    "--progress", 
    "--checksum")
  
  # Execute the command
  system(rclone_cmd, wait = TRUE)
}

# ---------------------------------------------------------
# Interpret '*' in dir name as not to be saved in Acacia
# ---------------------------------------------------------
identify_acacia_dirs = function(dir) {
  
  # Only output directories to be stored
  output_idx = grepl(dir$output, dir)
  dir$storage = dir[output_idx] %>%
    setdiff(dir$output) %>%
    unlist() %>%
    names()
  
  # Exluding any directories with '*' in the name
  #
  # NOTE: This should be any dirs with tons of files
  ignore_idx = grepl("\\*", dir[dir$storage])
  dir$ignore = dir$storage[ignore_idx]
  
  # Remove 'ignore' syntax from directory names
  dir = map(dir, ~str_remove(.x, "\\*"))
  
  return(dir)
}

# ---------------------------------------------------------
# Create a bucket on Acacia (if it doesn't already exist)
# ---------------------------------------------------------
create_acacia_link = function() {
  
  # Name of bucket we'll store to 
  bucket_name = get_bucket_name()
  
  # Vector of buckets that aleady exist
  buckets_exist = paste0("rclone lsf ", o$user, ":") %>%
    system(intern = TRUE) %>%
    substr(1, nchar(.) - 1)
  
  # Create the bucket if it doesn't already exist
  if (!bucket_name %in% buckets_exist)
    system(paste0("rclone mkdir ", get_bucket_path()))
}

# ---------------------------------------------------------
# Generate bucket name based on repo and user names
# ---------------------------------------------------------
get_bucket_name = function() {
  
  # Bucket name derived from repo name
  repo = basename(o$pth$code)
  
  # Combine version with user name, and remove invalid characters
  name = paste1(o$user, repo, o$output_version) %>%
    tolower() %>%
    # Only small letters, numbers, and dashes...
    str_replace_all("[^a-z0-9]", "-")
  
  return(name)
}

# ---------------------------------------------------------
# Concatenate full bucket path on Acacia
# ---------------------------------------------------------
get_bucket_path = function() {
  
  # Generate bucket name
  name = get_bucket_name()
  
  # Join with username for full path (the colon is important)
  path = paste(o$user, name, sep = ":")
  
  return(path)
}

