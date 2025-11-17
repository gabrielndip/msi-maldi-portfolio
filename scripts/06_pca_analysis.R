## PCA analysis and visualization script
##
## This script computes principal components analysis (PCA) on each of the
## processed datasets and visualises both the spatial **scores** and the
## feature **loadings**.
##
## - **Scores** are plotted as images to show the spatial distribution of the
##   principal components.
## - **Loadings** are plotted as spectra to show which m/z values contribute
##   most to each component, providing biological context.
##
## The script uses Cardinal's PCA implementation, which is built on
## the IRLBA algorithm for efficient computation.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Cardinal)
})
use_furrr <- FALSE
if (requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)) {
  use_furrr <- TRUE
}

### Configuration
# Set the number of principal components to calculate.
pca_ncomp <- 4

# The number of top loading values to label on the plots.
n_top_loadings <- 10
pca_figure_dir <- file.path("figures", "pca")
dir.create(pca_figure_dir, showWarnings = FALSE, recursive = TRUE)

# Determine which data list to use for analysis.
if (exists("msi_roi") && is.list(msi_roi) && length(msi_roi) > 0) {
  data_list <- msi_roi
  message("Using ROI-cropped datasets for PCA analysis.")
} else if (exists("msi_preprocessed") && is.list(msi_preprocessed) && length(msi_preprocessed) > 0) {
  data_list <- msi_preprocessed
  message("Using preprocessed datasets for PCA analysis.")
} else {
  stop("No suitable dataset list found. Please run preprocessing scripts first.")
}

# Create an empty list to hold PCA results.
msi_pca <- list()
msi_pca_mz <- list()
msi_pca_score_paths <- list()
msi_pca_loading_paths <- list()
msi_pca_loading_summary <- list()

render_pca_scores <- function(pca_result, dataset_name, outfile) {
  grDevices::png(outfile, width = 2200, height = 2200, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  nrows <- ceiling(pca_ncomp / 2)
  par(mfrow = c(nrows, 2), mar = c(2, 2, 3, 2), oma = c(0, 0, 3, 0))
  image(pca_result,
        superpose = FALSE,
        strip = FALSE,
        contrast.enhance = "histogram",
        normalize.image = "linear")
  mtext(paste0(dataset_name, ": PCA Scores"), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
}

render_pca_loadings <- function(pca_result, dataset_name, outfile) {
  grDevices::png(outfile, width = 2200, height = 2200, res = 220)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  nrows <- ceiling(pca_ncomp / 2)
  par(mfrow = c(nrows, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  loadings_df <- as.data.frame(loadings(pca_result))
  loadings_df$mz <- Cardinal::mz(pca_result)
  
  component_summary <- purrr::map_dfr(1:pca_ncomp, function(i) {
    pc_col <- paste0("PC", i)
    top_loadings <- loadings_df %>%
      arrange(desc(abs(.data[[pc_col]]))) %>%
      slice_head(n = n_top_loadings)
    
    plot(loadings_df$mz, loadings_df[[pc_col]], type = 'h',
         xlab = "m/z", ylab = "Loading",
         main = paste0("Loadings for ", pc_col))
    abline(h = 0, col = "red", lty = 2)
    
    text(top_loadings$mz, top_loadings[[pc_col]],
         labels = round(top_loadings$mz, 2),
         cex = 0.7, pos = 3, col = "blue")
    
    tibble(
      dataset = dataset_name,
      component = pc_col,
      top_mz = paste(round(head(top_loadings$mz, 5), 3), collapse = ", ")
    )
  })
  
  mtext(paste0(dataset_name, ": PCA Loadings"), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  component_summary
}

# Iterate over each dataset in the selected list.
mapper <- purrr::imap
cleanup_plan <- NULL
if (length(data_list) > 1 && use_furrr) {
  mapper <- furrr::future_imap
  cleanup_plan <- future::plan()
  future::plan(future::multisession, workers = min(4, length(data_list)))
}

mapper(data_list, function(msi_obj, name) {
  message(sprintf("Performing PCA on dataset '%s'...", name))
  
  # Perform PCA. Return a PCA2 object.
  pca_result <- tryCatch({
    PCA(msi_obj, ncomp = pca_ncomp)
  }, error = function(e) {
    stop(sprintf("Failed to compute PCA for dataset '%s': %s", name, e$message))
  })
  
  # Store and assign the result.
  msi_pca[[name]] <<- pca_result
  msi_pca_mz[[name]] <<- Cardinal::mz(msi_obj)
  assign(paste0(name, "_pca"), pca_result, envir = .GlobalEnv)
  
  #---------------------------------------------------------------------
  # 1. Plot PCA Scores
  #---------------------------------------------------------------------
  # These images show the spatial distribution of each principal component.
  tryCatch({
    score_file <- file.path(pca_figure_dir, paste0(name, "_scores.png"))
    render_pca_scores(pca_result, name, score_file)
    msi_pca_score_paths[[name]] <<- score_file
  }, error = function(e) {
    warning(sprintf("Could not plot PCA scores for dataset '%s': %s", name, e$message))
  })
  
  #---------------------------------------------------------------------
  # 2. Plot PCA Loadings
  #---------------------------------------------------------------------
  # These plots show which m/z values contribute most to each PC.
  # High positive or negative loadings indicate important features.
  tryCatch({
    loadings_file <- file.path(pca_figure_dir, paste0(name, "_loadings.png"))
    component_summary <- render_pca_loadings(pca_result, name, loadings_file)
    msi_pca_loading_paths[[name]] <<- loadings_file
    msi_pca_loading_summary[[name]] <<- component_summary
    
  }, error = function(e) {
    warning(sprintf("Could not plot PCA loadings for dataset '%s': %s", name, e$message))
  })
})

if (!is.null(cleanup_plan)) {
  future::plan(cleanup_plan)
}

# Assign the list of PCA results to the global environment.
assign("msi_pca", msi_pca, envir = .GlobalEnv)
assign("msi_pca_mz", msi_pca_mz, envir = .GlobalEnv)
assign("msi_pca_score_paths", msi_pca_score_paths, envir = .GlobalEnv)
assign("msi_pca_loading_paths", msi_pca_loading_paths, envir = .GlobalEnv)
assign(
  "msi_pca_loading_summary",
  bind_rows(msi_pca_loading_summary),
  envir = .GlobalEnv
)

message("PCA analysis complete. Results are in 'msi_pca' list and as individual variables.")

#-----------------------------------------------------------------------
# INTERPRETATION:
# The PCA score plots show the spatial distribution of variance. For example,
# PC1 might differentiate anatomical region A from B. To understand which
# molecules are responsible for this difference, we examine the loadings.
#
# In the loadings plot for PC1, m/z values with high positive loadings will
# be more abundant in the high-score regions of the PC1 image, while
# features with high negative loadings will be abundant in the low-score regions.
#
# NEXT STEP:
# Run 07_extra_analysis.R to perform more advanced analysis, or create a new
# script (e.g., 08_annotate_features.R) to match the important m/z values
# found here against a biological database.
#-----------------------------------------------------------------------
