###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

# # ---- R version check ----
# options(repos = c(CRAN = "https://cloud.r-project.org"))
# .libPaths(c("/software/projects/pawsey1163/johiolei/setonix/2025.08/r/4.4.1", .libPaths()))



# functions, helper functions and parameters 
source("R/sub_functions.R")
source("R/generic_function.R")
source("R/parameters_generic.R")

packages <- c("ggplot2", 
              "dplyr",
              "tibble",
              "tidyr", 
              "readr", 
              "purrr",       # uses pmap to loop through different scenarios
              "furrr",       # multisession i.e. distribute work across many cores
              "progressr",
              "patchwork",
              "sensitivity",
              "rsm",
              "randomForest",
              "ranger",
              "data.table",
              "rstudioapi")   # shows the progress


load_libraries(packages)