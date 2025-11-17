## 01_load_data.R
##
## Load MSI Datasets
##
## This script is the first step in the workflow. It locates all imzML
## files in the `data/` directory and loads them into R as `MSImageSet`
## objects using the `Cardinal::readMSIData()` function.
##
## The loaded datasets are placed into a list called `msi_data` and
## also into the global environment for easy access.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})
use_furrr <- FALSE
if (requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)) {
  use_furrr <- TRUE
}

## --------------------------------------------------------------------
## Configuration
## --------------------------------------------------------------------

# Define the directory containing the imzML/ibd pairs.
data_dir <- "data"

## --------------------------------------------------------------------
## Execution
## --------------------------------------------------------------------

message("--- Loading MSI datasets ---")

# Identify all .imzML files in the data directory.
imzml_files <- list.files(
  path = data_dir,
  pattern = "\\.imzML$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(imzml_files) == 0) {
  stop("No .imzML files were found in the 'data' directory.")
}

# Decide whether to load files in parallel. If there is more than one dataset
# and furrr/future are available, spawn a small multisession pool; otherwise
# fall back to sequential purrr::map to keep dependencies minimal.
mapper <- purrr::map
cleanup_plan <- NULL
if (length(imzml_files) > 1 && use_furrr) {
  mapper <- furrr::future_map
  cleanup_plan <- future::plan()
  future::plan(
    future::multisession,
    workers = min(4, length(imzml_files))
  )
}

# Load each imzML/ibd pair into an MSImageSet object.
# The .ibd file is automatically detected by Cardinal.
msi_data <- imzml_files %>%
  set_names(nm = ~ tools::file_path_sans_ext(basename(.))) %>%
  mapper(readMSIData)

if (!is.null(cleanup_plan)) {
  future::plan(cleanup_plan)
}

# Assign each data object into the global environment for convenience.
list2env(msi_data, envir = .GlobalEnv)

# Build a compact overview table so the Quarto report can render a
# dataset summary without re-reading the binary files.
msi_dataset_overview <- purrr::imap_dfr(msi_data, function(msi_obj, dataset_name) {
  mz_vals <- tryCatch(mz(msi_obj), error = function(e) NULL)
  tibble(
    dataset = dataset_name,
    pixels = ncol(spectra(msi_obj)),
    features = nrow(spectra(msi_obj)),
    mz_min = if (!is.null(mz_vals)) round(min(mz_vals, na.rm = TRUE), 4) else NA_real_,
    mz_max = if (!is.null(mz_vals)) round(max(mz_vals, na.rm = TRUE), 4) else NA_real_
  )
})
assign("msi_dataset_overview", msi_dataset_overview, envir = .GlobalEnv)

message("Successfully loaded datasets: ", paste(names(msi_data), collapse = ", "))

# Return the list invisibly.
invisible(msi_data)

## --------------------------------------------------------------------
# INTERPRETATION:
#
# This script populates your R environment with the raw MSI data. Each
# object (e.g., `Brain01_Bregma-1-46_centroid`) is a complex `MSImageSet`
# that contains all the spectra and spatial metadata for that sample.
#
# NEXT STEP:
# Run 02_preprocess_data.R to normalize the data and perform peak picking.
## --------------------------------------------------------------------
