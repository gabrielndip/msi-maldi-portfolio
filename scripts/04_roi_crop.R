## 04_roi_crop.R
##
## Crop to Region of Interest (ROI)
##
## This script provides an optional step to automatically crop the MSI
## datasets to the primary tissue region. This is useful for removing
## background pixels (e.g., the slide) to speed up downstream processing
## and focus the analysis on the tissue itself.
##
## The method is heuristic: it assumes that pixels with high Total Ion
## Current (TIC) correspond to tissue. It finds the bounding box around
## these high-TIC pixels and crops the data to that box.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})

## --------------------------------------------------------------------
## Configuration
## --------------------------------------------------------------------

# The quantile of pixel TICs to use for defining the tissue area.
# For example, 0.25 means the script will find the bounding box for all
# pixels with a TIC in the top 75% of all pixels.
# Adjust this value if your tissue occupies a smaller or larger portion
# of the total image area. A higher value (e.g., 0.5) will result in
# a tighter, more aggressive crop.
tic_quantile <- 0.25

## --------------------------------------------------------------------
## Setup
## --------------------------------------------------------------------
if (!exists("msi_preprocessed", inherits = TRUE)) {
  stop("msi_preprocessed not found. Please run 02_preprocess_data.R first.")
}

## --------------------------------------------------------------------
## ROI Cropping Function
## --------------------------------------------------------------------

#' Crop an MSI dataset to a bounding box defined by high-TIC pixels.
#'
#' @param proc_msi A processed MSImageSet object.
#' @param name The name of the dataset (for logging).
#' @param quantile The TIC quantile threshold for defining the ROI.
#' @return A cropped MSImageSet object, or NULL on failure.
crop_to_roi <- function(proc_msi, name, quantile = 0.25) {
  # Compute total ion current for each pixel.
  tic_values <- tryCatch(colSums(iData(proc_msi)), error = function(e) {
    warning(sprintf("Unable to compute TIC for %s: %s", name, e$message))
    return(NULL)
  })
  if (is.null(tic_values)) return(NULL)

  # Identify pixels with TIC above the specified quantile.
  threshold <- stats::quantile(tic_values, probs = quantile, na.rm = TRUE)
  high_tic_pixels <- which(tic_values > threshold)
  if (length(high_tic_pixels) == 0) {
    warning(sprintf("No pixels above TIC threshold for %s; skipping crop.", name))
    return(proc_msi) # Return original object if no crop is possible
  }

  # Get the coordinates of the high-TIC pixels and find their bounding box.
  pos <- pixelData(proc_msi)[high_tic_pixels, ]
  x_range <- range(pos$x, na.rm = TRUE)
  y_range <- range(pos$y, na.rm = TRUE)

  # Select all pixels within this bounding box.
  roi_pixels <- pixels(proc_msi, x >= x_range[1] & x <= x_range[2] &
                                 y >= y_range[1] & y <= y_range[2])
  
  # Subset the original object to these pixels.
  message(sprintf("  Cropping '%s': Original pixels = %d, ROI pixels = %d",
                  name, ncol(proc_msi), length(roi_pixels)))
  return(proc_msi[, roi_pixels])
}

## --------------------------------------------------------------------
## Execution
## --------------------------------------------------------------------

message("--- Cropping datasets to Region of Interest (ROI) ---")

# Apply the cropping function to each processed dataset.
msi_roi <- purrr::imap(msi_preprocessed, ~crop_to_roi(.x, .y, quantile = tic_quantile))

# Assign new objects to the global environment.
list2env(msi_roi, envir = .GlobalEnv)

message("ROI cropping complete. New objects are available (e.g., Brain01_..._proc_roi).")

invisible(msi_roi)

## --------------------------------------------------------------------
# INTERPRETATION:
#
# This script creates new `MSImageSet` objects, usually with a `_roi`
# suffix. These are smaller and contain only the main tissue region.
# Using these cropped objects can significantly speed up the subsequent
# computationally-intensive steps like segmentation and PCA.
#
# NEXT STEP:
# Run 05_ssc_segmentation.R to perform spatially-aware clustering on
# these focused, ROI-cropped datasets.
## --------------------------------------------------------------------
