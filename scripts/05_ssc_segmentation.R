## 05_segmentation.R
##
## Spatially-Aware and Non-Spatial Clustering
##
## This script performs segmentation to identify distinct molecular regions.
## It implements and compares three methods:
##
## 1.  **Spatial Shrunken Centroids (SSC)**: A sophisticated spatially-aware
##     method that performs feature selection while clustering.
## 2.  **Spatial K-Means (SKM)**: A faster spatially-aware method that
##     incorporates a spatial smoothing penalty into the k-means algorithm.
## 3.  **K-Means on PCA Scores**: A standard non-spatial method, used as a
##     baseline for comparison.
##
## The script first explores a range of `k` values for SSC, then generates
## a side-by-side comparison of the three methods for a fixed `k`.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})

## --------------------------------------------------------------------
## Configuration
## --------------------------------------------------------------------

# --- Spatial Method Parameters ---
# Neighborhood radius (r): Larger values increase spatial smoothing.
spatial_r <- 2

# Shrinkage parameter (s) for SSC: Higher values perform more aggressive
# feature selection.
ssc_s <- 15

# --- K-Value Parameters ---
# Range of k (segments) to test for the initial SSC exploration.
k_range <- 3:6

# The single 'k' to use for the final 3-method comparison.
comparison_k <- 4
segmentation_figure_dir <- file.path("figures", "segmentation")
dir.create(segmentation_figure_dir, showWarnings = FALSE, recursive = TRUE)

## --------------------------------------------------------------------
## Setup
## --------------------------------------------------------------------

# Choose the input datasets for segmentation.
if (exists("msi_roi") && is.list(msi_roi) && length(msi_roi) > 0) {
  message("Using ROI‑cropped datasets (msi_roi) for segmentation.")
  segmentation_input <- msi_roi
} else if (exists("msi_preprocessed") && is.list(msi_preprocessed) && length(msi_preprocessed) > 0) {
  message("Using preprocessed datasets (msi_preprocessed) for segmentation.")
  segmentation_input <- msi_preprocessed
} else {
  stop("No datasets available for segmentation. Run preprocessing scripts first.")
}

# Check if PCA results are available for k-means comparison.
if (!exists("msi_pca", inherits = TRUE) || !is.list(msi_pca)) {
  msi_pca <- list()
}

ensure_pca_baseline <- function(dataset_name, msi_obj, ncomp = max(comparison_k, 5)) {
  if (!dataset_name %in% names(msi_pca)) {
    message(sprintf("  Computing PCA for '%s' to support k-means baseline.", dataset_name))
    msi_pca[[dataset_name]] <<- tryCatch(
      PCA(msi_obj, ncomp = ncomp),
      error = function(e) {
        warning(sprintf("PCA failed for '%s': %s", dataset_name, e$message))
        NULL
      }
    )
  }
  msi_pca[[dataset_name]]
}

## --------------------------------------------------------------------
## 1. Explore k-values with Spatial Shrunken Centroids (SSC)
## --------------------------------------------------------------------

message(sprintf("Performing SSC segmentation for k = %s...", paste(k_range, collapse = ", ")))

msi_ssc <- list()
purrr::iwalk(segmentation_input, function(msi_obj, name) {
  msi_ssc[[name]] <<- list()
  for (k in k_range) {
    message(sprintf("  Applying SSC to '%s' (k = %d)...", name, k))
    res <- tryCatch({
      spatialShrunkenCentroids(msi_obj, method = "adaptive", r = spatial_r, s = ssc_s, k = k)
    }, error = function(e) {
      warning(sprintf("SSC failed for %s with k=%d: %s", name, k, e$message))
      return(NULL)
    })
    if (!is.null(res)) msi_ssc[[name]][[paste0("k", k)]] <<- res
  }
})

# Plot the SSC results for each k to allow visual comparison and capture PNGs.
message("Plotting SSC results for different k values...")
msi_segmentation_plots <- list()
msi_segmentation_overview <- list()

save_ssc_grid <- function(results_list, dataset_name) {
  outfile <- file.path(segmentation_figure_dir, paste0(dataset_name, "_ssc_grid.png"))
  grDevices::png(outfile, width = 2200, height = 800, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, length(results_list)), oma = c(0, 0, 3, 0))
  purrr::iwalk(results_list, ~plot(.x, main = .y))
  mtext(paste0(dataset_name, ": SSC Segmentation Results (Varying k)"), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  outfile
}

purrr::iwalk(msi_ssc, function(results_list, name) {
  if (length(results_list) == 0) return()
  grid_path <- save_ssc_grid(results_list, name)
  msi_segmentation_plots[[paste0(name, "_ssc_grid")]] <<- grid_path
  msi_segmentation_overview[[name]] <<- tibble(
    dataset = name,
    k_sweep = paste(names(results_list), collapse = ", ")
  )
})

## --------------------------------------------------------------------
## 2. Compare Segmentation Methods (SSC vs. SKM vs. K-Means)
## --------------------------------------------------------------------

message(sprintf("\nPerforming 3-way method comparison for k=%d...", comparison_k))
msi_skm <- list() # To store SKM results

purrr::iwalk(segmentation_input, function(msi_obj, name) {
  
  # --- Data Check ---
  ssc_result <- msi_ssc[[name]][[paste0("k", comparison_k)]]
  pca_result <- ensure_pca_baseline(name, msi_obj)
  if (is.null(ssc_result) || is.null(pca_result)) {
    warning(sprintf("Skipping comparison for '%s' due to missing SSC (for k=%d) or PCA results.", name, comparison_k))
    return()
  }
  message(sprintf("  Comparing methods for '%s'...", name))
  
  # --- Run Additional Models ---
  # a) Spatial K-Means (SKM)
  skm_result <- tryCatch({
    spatialKMeans(msi_obj, r = spatial_r, k = comparison_k, method = "adaptive")
  }, error = function(e) {
    warning(sprintf("SKM failed for %s: %s", name, e$message)); return(NULL)
  })
  if (!is.null(skm_result)) msi_skm[[name]] <<- skm_result
  
  # b) Non-spatial K-Means on PCA scores
  set.seed(123)
  kmeans_res <- kmeans(scores(pca_result), centers = comparison_k)
  pData(msi_obj)$kmeans_cluster <- factor(kmeans_res$cluster)
  
  # --- Plot 3-Way Comparison ---
  comparison_file <- file.path(segmentation_figure_dir, paste0(name, "_method_comparison.png"))
  grDevices::png(comparison_file, width = 2200, height = 800, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 3), oma = c(0, 0, 3, 0))
  
  # Plot 1: SSC Result
  plot(ssc_result, main = "Method: SSC")
  
  # Plot 2: SKM Result
  if (!is.null(skm_result)) {
    plot(skm_result, main = "Method: SKM")
  } else {
    plot.new(); title("SKM Failed")
  }
  
  # Plot 3: K-means Result
  plot(msi_obj, .color = ~ kmeans_cluster, main = "Method: K-Means (non-spatial)")
  
  mtext(paste0(name, ": Segmentation Method Comparison (k=", comparison_k, ")"), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  msi_segmentation_plots[[paste0(name, "_comparison")]] <<- comparison_file
  existing_summary <- msi_segmentation_overview[[name]]
  comparison_row <- tibble(
    dataset = name,
    k_sweep = if (!is.null(existing_summary$k_sweep)) existing_summary$k_sweep else paste(names(msi_ssc[[name]]), collapse = ", "),
    comparison_k = comparison_k,
    has_skm = !is.null(skm_result)
  )
  msi_segmentation_overview[[name]] <<- comparison_row
})

# --- Save results to environment ---
assign("msi_ssc", msi_ssc, envir = .GlobalEnv)
assign("msi_skm", msi_skm, envir = .GlobalEnv)
assign("msi_pca", msi_pca, envir = .GlobalEnv)
assign("msi_segmentation_plots", msi_segmentation_plots, envir = .GlobalEnv)
assign(
  "msi_segmentation_summary",
  bind_rows(msi_segmentation_overview),
  envir = .GlobalEnv
)
message("\nSegmentation analysis complete.")

## --------------------------------------------------------------------
# INTERPRETATION:
#
# - The first set of plots (SSC with varying k) helps you choose a `k` that
#   best matches the biological structures in your tissue.
#
# - The 3-way comparison plot demonstrates the impact of spatial information.
#   - **K-Means** is often noisy (a "salt-and-pepper" look) because it
#     treats every pixel independently.
#   - **SKM** and **SSC** produce smoother, more contiguous regions because
#     they consider neighboring pixels. SSC is often even more refined as
#     it also performs feature selection simultaneously.
#
# NEXT STEP:
# Run scripts 06, 07, and 08 to identify the molecules that define these
# segments and attempt to annotate them.
## --------------------------------------------------------------------
