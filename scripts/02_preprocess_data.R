## 02_preprocess_data.R
##
## Preprocess MSI Datasets
##
## This script applies a standard preprocessing pipeline to each of the
## raw `MSImageSet` objects loaded by the previous script. This is a
## critical step to clean the data and make spectra comparable.
##
## The pipeline consists of:
## 1.  **Normalization**: To correct for technical variations in signal
##     intensity across the tissue.
## 2.  **Peak Picking**: To identify m/z values that represent real signals.
## 3.  **Peak Alignment**: To correct for small mass shifts between pixels.
## 4.  **Peak Filtering**: To remove peaks that are very rare.
##
## The processed datasets are saved in a list called `msi_preprocessed`.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})

## --------------------------------------------------------------------
## Setup
## --------------------------------------------------------------------
if (!exists("msi_data", inherits = TRUE)) {
  stop("msi_data not found. Please run 01_load_data.R first.")
}

## --------------------------------------------------------------------
## Preprocessing Function
## --------------------------------------------------------------------

#' Apply a standard preprocessing pipeline to an MSImageSet.
#'
#' @param msi An MSImageSet object.
#' @return A processed MSImageSet object.
preprocess_msi <- function(msi) {
  msi %>%
    # Normalize by Total Ion Current (TIC). This corrects for differences
    # in the total signal intensity between pixels.
    normalize(method = "tic") %>%
    
    # Pick peaks using the Median Absolute Deviation (MAD) method. This
    # identifies peaks that stand out from the spectral noise. SNR is the
    # signal-to-noise ratio threshold.
    peakPick(method = "mad", SNR = 5) %>%
    
    # Align peaks across all pixels to a common m/z axis. The tolerance
    # defines how far a peak can be to be considered the "same" peak.
    peakAlign(tolerance = 0.5, units = "mz") %>%
    
    # Filter peaks. Here, we remove peaks that are present in fewer than
    # 1% of the pixels (freq.min = 0.01).
    peakFilter(freq.min = 0.01, rm.zero = TRUE) %>%
    
    # Execute the queued processing steps.
    process()
}

## --------------------------------------------------------------------
## Execution
## --------------------------------------------------------------------

message("--- Preprocessing datasets ---")

# Apply the preprocessing pipeline to each dataset.
msi_preprocessed <- msi_data %>%
  purrr::imap(~{
    message("Processing: ", .y)
    preprocess_msi(.x)
  }) %>%
  # Append _proc suffix to names to distinguish from raw data.
  set_names(~ paste0(names(msi_data), "_proc"))

# Assign processed datasets into the global environment.
list2env(msi_preprocessed, envir = .GlobalEnv)

message(
  "Preprocessing complete. New objects: ",
  paste(names(msi_preprocessed), collapse = ", ")
)

invisible(msi_preprocessed)

## --------------------------------------------------------------------
# INTERPRETATION:
#
# This script transforms the raw data into a processed form suitable for
# analysis. The key change is that the data is no longer continuous spectra,
# but a set of discrete, aligned peaks. This reduces noise and data size,
# and makes all the downstream analyses (PCA, clustering) possible.
#
# NEXT STEP:
# Run 03_explore_data.R to perform quality control on the processed data
# and to see the first images of biologically-relevant molecules.
## --------------------------------------------------------------------
