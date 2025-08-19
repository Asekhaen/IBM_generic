###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

# ---- R version check ----

# R versions for which this project has been tested and is stable
stable_versions = c("4.4.0", "4.4.1")

# R versions for which this project is stable (as a string)
stable_str = paste(stable_versions, collapse = ", ")

# Get details of R version currently running
version_info = R.Version()

# Construct version number from list details
version_num = paste0(version_info$major, ".",  version_info$minor)

# Throw an error if this R version is unsuitable
if (!version_num %in% stable_versions)
  stop("This software is stable with R version(s): ", stable_str,
       " (currently running ", version_num, ")")

# Clear global environment
rm(list = ls())

# ---- Source files ----

# Scripts that should not be sourced
no_src = c(
  "launch.R", 
  "dependencies.R", 
  "submit.R") 

# All R files, and those to source
all_files = list.files(pattern = ".+\\.R$")
src_files = setdiff(all_files, no_src)

# Source each of these files
for (file in src_files)
  source(file)

# ---- Define packages ----

# Complete list of all packages required for this project
packages = c(
  # "tidyverse",     # Includes ggplot2, dplyr, tidyr (www.tidyverse.org/packages/)
  "ggplot2",       # TIDYVERSE
  "dplyr",         # TIDYVERSE
  "tibble",        # TIDYVERSE
  "data.table",    # Next generation dataframes
  "dtplyr",        # Syntax of dplyr with datatable speed
  "tidyr",         # TIDYVERSE
  "readr",         # TIDYVERSE
  "purrr",         # TIDYVERSE
  "stringr",       # TIDYVERSE
  "forcats",       # TIDYVERSE
  "lubridate",     # TIDYVERSE
  "magrittr",      # Additional pipe operators, such as %<>%
  "yaml",          # Read YAML markup files
  "wrapr",         # Convenience functions, such as qc
  "parallel",      # Multicore version of lapply
  "rlist",         # List-related helper functions, such as list.remove
  "purrr",         # TODO: Write description
  "truncnorm",     # TODO: Write description
  "parallel",      # Multicore version of lapply
  "pals",          # Colour palettes
  "progress")      # Stylish progress bars

# ---- Install and/or load packages with pacman ----

message("* Installing required packages")

# Check whether pacman itself has been installed
pacman_installed = "pacman" %in% rownames(installed.packages())

# If not, install it
if (!pacman_installed) 
  install.packages("pacman")

# Load pacman
library(pacman) 

# Load all required packages, installing them if required
pacman::p_load(char = packages)

# ---- Redefine or unmask particular functions ----

# Unmask certain functions otherwise overwritten
select  = dplyr::select
filter  = dplyr::filter
rename  = dplyr::rename
recode  = dplyr::recode
count   = dplyr::count
predict = stats::predict

# ---- Tidy up ----

# Tidy up
if (interactive()) clc()  # Clear console
if (interactive()) clf()  # Close figures

