



# Package
library(terra)
library(glue)
library(terra)
library(hdar)      # to access Copernicus data using wekeo services
library(jsonlite)


root_dir   <- getwd()
input_case <- file.path(root_dir, "input")
functions_case <- file.path(root_dir, "functions")
output_case <- file.path(root_dir, "output")

dir.create(input_case, recursive = TRUE, showWarnings = FALSE)
dir.create(functions_case, recursive = TRUE, showWarnings = FALSE)
dir.create(output_case, recursive = TRUE, showWarnings = FALSE)



function_files <- list.files(functions_case, pattern = "\\.R$", full.names = TRUE)
