## 03_explore_data.R
##
## Quality Control and Exploratory Data Analysis
##
## This script produces a set of diagnostic visualizations for each raw
## and processed dataset. The goal is to assess data quality and get a
## first look at the most prominent spatial patterns.
##
## The script generates three key sets of plots:
## 1.  **Total Ion Current (TIC) Images**: To check for overall signal
##     intensity and identify the tissue area.
## 2.  **Spatially Ranked Ion Images**: Instead of picking an arbitrary ion,
##     this script now identifies the top 3 features with the highest
##     spatial structure (using Moran's I statistic) and displays their
##     ion images. This immediately highlights biologically relevant patterns.
## 3.  **Mass Spectra**: Shows the m/z profile for a randomly selected pixel.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})

## --------------------------------------------------------------------
## Configuration
## --------------------------------------------------------------------
set.seed(123) # For reproducible pixel selection
N_SPATIAL_FEATURES <- 3 # Number of top spatial features to plot

## --------------------------------------------------------------------
## Setup
## --------------------------------------------------------------------
if (!exists("msi_data", inherits = TRUE) || !exists("msi_preprocessed", inherits = TRUE)) {
  stop("Raw 'msi_data' or 'msi_preprocessed' not found. Run 01 and 02 scripts first.")
}

## --------------------------------------------------------------------
## Plotting Functions
## --------------------------------------------------------------------

#' Plot TIC images for a raw and processed dataset.
plot_tic_images <- function(raw_msi, proc_msi) {
  plot(raw_msi, main = "Raw TIC")
  plot(proc_msi, main = "Processed TIC")
}

#' Plot spectra for a single random pixel.
plot_spectra <- function(raw_msi, proc_msi) {
  pixel_count <- ncol(iData(raw_msi))
  pixel_idx <- if (is.na(pixel_count) || pixel_count < 1) 1 else sample(pixel_count, 1)
  
  plot(raw_msi, pixel = pixel_idx, main = paste("Raw Spectrum - Pixel", pixel_idx))
  plot(proc_msi, pixel = pixel_idx, main = paste("Processed Spectrum - Pixel", pixel_idx))
}

#' Find and plot the most spatially-structured ion images.
#'
#' @param msi_obj The processed MSImageSet.
#' @param n The number of top features to plot.
plot_spatially_ranked_images <- function(msi_obj, n = 3) {
  message(sprintf("  Finding top %d spatially-ranked features...", n))
  
  # Use spatialFeatures with local Moran's I ('rI') to rank features.
  # This is computationally intensive. We use a subset of pixels for speed.
  top_features <- tryCatch({
    spatialFeatures(msi_obj, method = "rI", r = 2, top = n)
  }, error = function(e) {
    warning("spatialFeatures failed. Falling back to top intensity features. Error: ", e$message)
    # Fallback: if spatialFeatures fails, use top N by total intensity.
    top_intensity_idx <- order(rowSums(iData(msi_obj)), decreasing = TRUE)[1:n]
    features(msi_obj, mz = mz(msi_obj)[top_intensity_idx])
  })
  
  if (length(top_features) == 0) {
    warning("Could not identify top features to plot.")
    replicate(n, plot.new()) # Plot empty frames to keep layout consistent
    return()
  }
  
  # Plot the ion image for each top feature
  for (i in 1:length(top_features)) {
    f <- top_features[i]
    mz_val <- mz(f)
    image(f, main = sprintf("Top Spatial Feature %d\nm/z = %.2f", i, mz_val))
  }
  
  # If we plotted fewer than n features, fill the rest with empty plots
  if (length(top_features) < n) {
    replicate(n - length(top_features), plot.new())
  }
}

## --------------------------------------------------------------------
## Main Execution Loop
## --------------------------------------------------------------------

# Iterate over each dataset to generate and save plots.
purr::iwalk(msi_data, function(raw_msi, name) {
  proc_name <- paste0(name, "_proc")
  if (!proc_name %in% names(msi_preprocessed)) {
    warning("No processed dataset found for ", name, " – skipping.")
    return()
  }
  proc_msi <- msi_preprocessed[[proc_name]]
  
  message(sprintf("\n--- Generating QC plots for: %s ---", name))
  
  # Set up plot layout: 2 rows, 3 columns
  # Row 1: TIC images and a spectrum
  # Row 2: Top 3 spatially-ranked ion images
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 5, 6, 6, 7, 8), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 3, 2), oma = c(0, 0, 3, 0))
  
  # --- Generate Plots ---
  tryCatch({
    # Plot TIC (Raw vs Processed)
    plot(raw_msi, main = "Raw TIC"); plot(proc_msi, main = "Processed TIC")
    
    # Plot Spectrum (Raw vs Processed)
    pixel_idx <- if (ncol(iData(raw_msi)) < 1) 1 else sample(ncol(iData(raw_msi)), 1)
    plot(raw_msi, pixel = pixel_idx, main = paste("Raw Spectrum", pixel_idx))
    plot(proc_msi, pixel = pixel_idx, main = paste("Processed Spectrum", pixel_idx))
    
    # Plot Spatially Ranked Features from the processed data
    plot_spatially_ranked_images(proc_msi, n = N_SPATIAL_FEATURES)
    
    # Add an overall title
    mtext(paste("QC Report:", name), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
    
  }, error = function(e) {
    warning(sprintf("Failed to generate full QC plot for %s: %s", name, e$message))
  })
})

message("\nExploratory analysis complete.")

## --------------------------------------------------------------------
# INTERPRETATION:
#
# - The TIC images give a general overview of the data quality. The
#   processed TIC should look smoother than the raw.
#
# - The "Top Spatial Feature" plots are the most important. They show
#   the ion images for m/z values that have strong, non-random spatial
#   distributions. These are excellent candidates for being biomarkers
#   that define specific anatomical regions.
#
# NEXT STEP:
# Run 04_roi_crop.R to crop the data to the tissue area, followed by
# 05_ssc_segmentation.R to formally segment the regions hinted at here.
## --------------------------------------------------------------------
