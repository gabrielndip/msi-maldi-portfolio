## 03_explore_data.R
##
## Quality Control and Exploratory Data Analysis
##
## This script produces and saves a set of diagnostic visualizations for each
## raw and processed dataset. The goal is to assess data quality and get a
## first look at the most prominent spatial patterns.
##
## For each dataset, the script generates three separate plot panels:
## 1.  **TIC Panel**: TIC images for raw vs. processed data.
## 2.  **Spectra Panel**: Mass spectra for a random pixel for raw vs. processed data.
## 3.  **Spatial Features Panel**: Ion images for the top 3 features with the
##     highest spatial structure (using Moran's I).
##
## These plots are saved to the 'figures/qc' directory and are intended to be
## displayed in the final knitted HTML report.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})

## --------------------------------------------------------------------
## Configuration
## --------------------------------------------------------------------
set.seed(123) # For reproducible pixel selection
N_SPATIAL_FEATURES <- 3 # Number of top spatial features to plot
qc_figure_dir <- file.path("figures", "qc")
dir.create(qc_figure_dir, showWarnings = FALSE, recursive = TRUE)

## --------------------------------------------------------------------
## Setup
## --------------------------------------------------------------------
if (!exists("msi_data", inherits = TRUE) || !exists("msi_preprocessed", inherits = TRUE)) {
  stop("Raw 'msi_data' or 'msi_preprocessed' not found. Run 01 and 02 scripts first.")
}

## --------------------------------------------------------------------
## Plotting and Saving Functions
## --------------------------------------------------------------------

#' Save a panel comparing Raw and Processed TIC images.
save_tic_panel <- function(raw_msi, proc_msi, dataset_name) {
  outfile <- file.path(qc_figure_dir, paste0(dataset_name, "_qc_tic.png"))
  grDevices::png(outfile, width = 1200, height = 600, res = 150)
  on.exit(grDevices::dev.off())

  par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))
  
  plot(raw_msi, main = paste(dataset_name, "- Raw TIC"))
  
  if (!is.null(proc_msi) && nrow(proc_msi) > 0) {
    plot(proc_msi, main = "Processed TIC")
  } else {
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    title(main = "Processed TIC")
    text(1, 1, "No features remaining after preprocessing.", cex = 1.2)
  }
  
  return(outfile)
}

#' Save a panel comparing Raw and Processed spectra for a random pixel.
save_spectra_panel <- function(raw_msi, proc_msi, dataset_name) {
  outfile <- file.path(qc_figure_dir, paste0(dataset_name, "_qc_spectra.png"))
  grDevices::png(outfile, width = 1200, height = 600, res = 150)
  on.exit(grDevices::dev.off())
  
  pixel_idx <- if (ncol(spectra(raw_msi)) < 1) 1 else sample(ncol(spectra(raw_msi)), 1)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))
  
  plot(raw_msi, i = pixel_idx, main = paste(dataset_name, "- Raw Spectrum"))
  
  if (!is.null(proc_msi) && nrow(proc_msi) > 0) {
    plot(proc_msi, i = pixel_idx, main = "Processed Spectrum")
  } else {
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    title(main = "Processed Spectrum")
    text(1, 1, "No features remaining after preprocessing.", cex = 1.2)
  }

  return(outfile)
}

#' Find, plot, and save the most spatially-structured ion images.
save_spatial_panel <- function(msi_obj, dataset_name, n = 3) {
  outfile <- file.path(qc_figure_dir, paste0(dataset_name, "_qc_spatial.png"))
  grDevices::png(outfile, width = 1800, height = 600, res = 150)
  on.exit(grDevices::dev.off())

  if (is.null(msi_obj) || nrow(msi_obj) == 0) {
    warning("Skipping spatial analysis for ", dataset_name, ": No features in processed data.", call. = FALSE)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    title(main = paste(dataset_name, "- Top Spatial Features"))
    text(1, 1, "No features available for spatial analysis.", cex = 1.5)
    return(list(path = outfile, spatial_mz = NULL))
  }

  message(sprintf("  Finding top %d spatially-ranked features for %s...", n, dataset_name))

  # spatialStats was removed in newer Cardinal releases. Prefer it when available,
  # otherwise fall back to simple intensity-based ranking.
  spatial_stats_fn <- get0("spatialStats", envir = asNamespace("Cardinal"), inherits = FALSE)

  feature_idx <- integer(0)

  if (is.function(spatial_stats_fn)) {
    spatial_stats_result <- try(
      spatial_stats_fn(msi_obj, method = "moran", type = "local", r = 2),
      silent = TRUE
    )

    if (!inherits(spatial_stats_result, "try-error")) {
      spatial_stats_df <- as.data.frame(spatial_stats_result) %>%
        arrange(desc(MoranI)) %>%
        head(n)
      feature_idx <- match(spatial_stats_df$mz, mz(msi_obj))
    } else {
      warning("spatialStats failed for ", dataset_name, ". Falling back to top intensity.", call. = FALSE)
    }
  } else {
    warning("Cardinal::spatialStats not available; using intensity ranking instead.", call. = FALSE)
  }

  if (length(feature_idx) == 0 || all(is.na(feature_idx))) {
    feature_idx <- order(rowSums(spectra(msi_obj)), decreasing = TRUE)
  }

  feature_idx <- unique(stats::na.omit(feature_idx))
  feature_idx <- head(feature_idx, n)

  if (length(feature_idx) == 0) {
    warning("Could not identify top features to plot for ", dataset_name, call. = FALSE)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    title(main = paste(dataset_name, "- Top Spatial Features"))
    text(1, 1, "Could not identify any top spatial features.", cex = 1.5)
    return(list(path = outfile, spatial_mz = NULL))
  }

  par(mfrow = c(1, length(feature_idx)), mar = c(2, 2, 5, 1))

  mz_values <- mz(msi_obj)[feature_idx]

  purrr::iwalk(feature_idx, function(idx, panel_idx) {
    image(
      msi_obj,
      i = idx,
      main = sprintf("Top Feature %d\nm/z = %.2f", panel_idx, mz_values[panel_idx])
    )
  })

  return(list(path = outfile, spatial_mz = unique(round(mz_values, 4))))
}


## --------------------------------------------------------------------
## Main Execution Loop
## --------------------------------------------------------------------

msi_qc_plots <- list()
msi_qc_spatial_features <- list()

for (name in names(msi_data)) {
  raw_msi <- msi_data[[name]]
  proc_name <- paste0(name, "_proc")
  
  if (!proc_name %in% names(msi_preprocessed)) {
    warning("No processed dataset found for ", name, " – skipping QC.")
    next
  }
  proc_msi <- msi_preprocessed[[proc_name]]
  
  message(sprintf("\n--- Generating QC plots for: %s ---", name))
  
  tryCatch({
    tic_path <- save_tic_panel(raw_msi, proc_msi, name)
    spectra_path <- save_spectra_panel(raw_msi, proc_msi, name)
    spatial_info <- save_spatial_panel(proc_msi, name, n = N_SPATIAL_FEATURES)
    
    msi_qc_plots[[name]] <- list(
      tic = tic_path,
      spectra = spectra_path,
      spatial = spatial_info$path
    )
    msi_qc_spatial_features[[name]] <- spatial_info$spatial_mz
    
  }, error = function(e) {
    warning(sprintf("Failed to generate QC plots for %s: %s", name, e$message))
  })
}

message("\nExploratory analysis complete. QC plots saved to 'figures/qc/'.")
assign("msi_qc_plots", msi_qc_plots, envir = .GlobalEnv)
assign("msi_qc_spatial_features", msi_qc_spatial_features, envir = .GlobalEnv)


## --------------------------------------------------------------------
# INTERPRETATION:
#
# - The TIC images give a general overview of the data quality.
# - The "Top Spatial Feature" plots are the most important, showing
#   biomarkers that define specific anatomical regions.
#
# NEXT STEP:
# The Quarto document will now display these generated plots. Following this,
# run 04_roi_crop.R to crop the data to the tissue area.
## --------------------------------------------------------------------
