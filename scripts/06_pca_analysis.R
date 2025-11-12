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

### Configuration
# Set the number of principal components to calculate.
pca_ncomp <- 4

# The number of top loading values to label on the plots.
n_top_loadings <- 10

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

# Iterate over each dataset in the selected list.
purrr::imap(data_list, function(msi_obj, name) {
  message(sprintf("Performing PCA on dataset '%s'...", name))
  
  # Perform PCA. Return a PCA2 object.
  pca_result <- tryCatch({
    PCA(msi_obj, ncomp = pca_ncomp)
  }, error = function(e) {
    stop(sprintf("Failed to compute PCA for dataset '%s': %s", name, e$message))
  })
  
  # Store and assign the result.
  msi_pca[[name]] <<- pca_result
  assign(paste0(name, "_pca"), pca_result, envir = .GlobalEnv)
  
  #---------------------------------------------------------------------
  # 1. Plot PCA Scores
  #---------------------------------------------------------------------
  # These images show the spatial distribution of each principal component.
  tryCatch({
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    nrows <- ceiling(pca_ncomp / 2)
    par(mfrow = c(nrows, 2), mar = c(2, 2, 3, 2), oma = c(0, 0, 3, 0))
    image(pca_result,
          superpose = FALSE,
          strip = FALSE,
          contrast.enhance = "histogram",
          normalize.image = "linear")
    mtext(paste0(name, ": PCA Scores"), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
  }, error = function(e) {
    warning(sprintf("Could not plot PCA scores for dataset '%s': %s", name, e$message))
  })
  
  #---------------------------------------------------------------------
  # 2. Plot PCA Loadings
  #---------------------------------------------------------------------
  # These plots show which m/z values contribute most to each PC.
  # High positive or negative loadings indicate important features.
  tryCatch({
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    nrows <- ceiling(pca_ncomp / 2)
    par(mfrow = c(nrows, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
    
    loadings_df <- as.data.frame(loadings(pca_result))
    loadings_df$mz <- Cardinal::mz(pca_result)
    
    for (i in 1:pca_ncomp) {
      pc_col <- paste0("PC", i)
      
      # Get top n loadings by absolute value for labeling
      top_loadings <- loadings_df %>%
        arrange(desc(abs(.data[[pc_col]]))) %>%
        slice_head(n = n_top_loadings)
      
      plot(loadings_df$mz, loadings_df[[pc_col]], type = 'h',
           xlab = "m/z", ylab = "Loading",
           main = paste0("Loadings for ", pc_col))
      abline(h = 0, col = "red", lty = 2)
      
      # Add labels for top loadings
      text(top_loadings$mz, top_loadings[[pc_col]],
           labels = round(top_loadings$mz, 2),
           cex = 0.7, pos = 3, col = "blue")
    }
    
    mtext(paste0(name, ": PCA Loadings"), side = 3, line = 1, outer = TRUE, font = 2, cex = 1.2)
    
  }, error = function(e) {
    warning(sprintf("Could not plot PCA loadings for dataset '%s': %s", name, e$message))
  })
})

# Assign the list of PCA results to the global environment.
assign("msi_pca", msi_pca, envir = .GlobalEnv)

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
