## Extra analyses script
##
## This script demonstrates several additional analyses on the processed
## MSI datasets, focusing on biological pattern discovery. It covers:
##
## 1.  **Manifold Learning (t-SNE & UMAP)**: To visualize the complex,
##     non-linear relationships between pixel spectra.
## 2.  **Supervised Classification**: To identify which m/z features are most
##     predictive of the spatial regions found during segmentation.
## 3.  **Co-localization Analysis**: To build networks of ions that are
##     frequently found together.

suppressPackageStartupMessages({
  library(Cardinal)
  library(tidyverse)
})

## --------------------------------------------------------------------
## Configuration & Setup
## --------------------------------------------------------------------

# Load optional packages, warning if they are not available.
packages <- c("Rtsne", "uwot", "randomForest", "igraph")
loaded_pkgs <- sapply(packages, function(pkg) {
  suppressPackageStartupMessages({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("Package '%s' not installed; skipping related analysis.", pkg))
      return(FALSE)
    }
    return(TRUE)
  })
})

# Determine which datasets to analyse.
if (exists("msi_roi") && is.list(msi_roi) && length(msi_roi) > 0) {
  data_list <- msi_roi
  message("Using ROI-cropped datasets for extra analyses.")
} else if (exists("msi_preprocessed") && is.list(msi_preprocessed) && length(msi_preprocessed) > 0) {
  data_list <- msi_preprocessed
  message("Using preprocessed datasets for extra analyses.")
} else {
  stop("No suitable dataset list found. Please run preprocessing scripts first.")
}

# Check for segmentation results, which will be used as labels.
seg_available <- exists("msi_ssc") && is.list(msi_ssc) && length(msi_ssc) > 0

## --------------------------------------------------------------------
## Analysis Functions
## --------------------------------------------------------------------

#' Extract top variable features and labels for analysis.
#'
#' @param msi_obj The MSImageSet to process.
#' @param dataset_name The name of the dataset (for logging).
#' @param nfeat The number of top features to select based on variance.
#' @return A list containing the feature matrix (pixels x features) and labels.
extract_features <- function(msi_obj, dataset_name, nfeat = 200) {
  # Extract intensity matrix (features × pixels)
  int_mat <- tryCatch(as.matrix(iData(msi_obj)), error = function(e) {
    warning(sprintf("Failed to extract intensity matrix for '%s': %s", dataset_name, e$message))
    return(NULL)
  })
  if (is.null(int_mat)) return(NULL)
  
  # Transpose to get pixels × features
  mat <- t(int_mat)
  
  # Filter to top 'nfeat' most variable features
  var_feat <- apply(mat, 2, var)
  nf <- min(nfeat, length(var_feat))
  top_idx <- order(var_feat, decreasing = TRUE)[seq_len(nf)]
  mat_filtered <- mat[, top_idx, drop = FALSE]
  
  # Retrieve segmentation labels if available
  labels <- NULL
  if (seg_available && !is.null(msi_ssc[[dataset_name]])) {
    labels <- tryCatch(modelData(msi_ssc[[dataset_name]])$class, error = function(e) {
      warning(sprintf("Could not extract class labels for '%s': %s", dataset_name, e$message))
      return(NULL)
    })
  }
  return(list(matrix = mat_filtered, labels = labels))
}


#' Perform and plot t-SNE and UMAP.
#'
#' @param features A matrix of features (pixels x features).
#' @param labels A vector of class labels for coloring points.
#' @param dataset_name The name of the dataset for plot titles.
perform_manifold <- function(features, labels, dataset_name) {
  # t-SNE analysis
  if (loaded_pkgs["Rtsne"]) {
    tryCatch({
      set.seed(123)
      tsne_res <- Rtsne::Rtsne(features, dims = 2, perplexity = 30, verbose = FALSE, check_duplicates = FALSE)
      tsne_df <- data.frame(tsne1 = tsne_res$Y[,1], tsne2 = tsne_res$Y[,2], class = factor(labels))
      p_tsne <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = class)) +
        geom_point(size = 1, alpha = 0.7) +
        scale_color_viridis_d(option = "plasma", na.value = "grey50") +
        labs(title = paste0(dataset_name, ": t-SNE"), color = "Class") +
        theme_minimal()
      print(p_tsne)
    }, error = function(e) warning(sprintf("t-SNE failed for '%s': %s", dataset_name, e$message)))
  }
  # UMAP analysis
  if (loaded_pkgs["uwot"]) {
    tryCatch({
      set.seed(123)
      umap_res <- uwot::umap(features, n_components = 2)
      umap_df <- data.frame(UMAP1 = umap_res[,1], UMAP2 = umap_res[,2], class = factor(labels))
      p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = class)) +
        geom_point(size = 1, alpha = 0.7) +
        scale_color_viridis_d(option = "plasma", na.value = "grey50") +
        labs(title = paste0(dataset_name, ": UMAP"), color = "Class") +
        theme_minimal()
      print(p_umap)
    }, error = function(e) warning(sprintf("UMAP failed for '%s': %s", dataset_name, e$message)))
  }
}


#' Train a Random Forest and plot feature importance.
#'
#' @param features A matrix of features (pixels x features).
#' @param labels A vector of class labels.
#' @param dataset_name The name of the dataset.
#' @return The trained randomForest model.
perform_classification <- function(features, labels, dataset_name) {
  if (!loaded_pkgs["randomForest"]) return(NULL)
  if (is.null(labels)) {
    warning(sprintf("No labels available for supervised classification on '%s'.", dataset_name))
    return(NULL)
  }
  tryCatch({
    labels <- factor(labels)
    set.seed(123)
    train_idx <- sample(seq_len(nrow(features)), size = floor(0.8 * nrow(features)))
    
    rf_model <- randomForest::randomForest(x = features[train_idx,], y = labels[train_idx],
                                           xtest = features[-train_idx,], ytest = labels[-train_idx],
                                           ntree = 100)
    
    acc <- mean(rf_model$test$predicted == rf_model$test$y)
    message(sprintf("Random Forest accuracy on test set for '%s': %.2f%%", dataset_name, 100 * acc))
    
    # --- Plot Feature Importance ---
    imp_df <- randomForest::importance(rf_model) %>%
      as.data.frame() %>%
      rownames_to_column("mz") %>%
      arrange(desc(MeanDecreaseGini)) %>%
      slice_head(n = 20) %>%
      mutate(mz = fct_reorder(mz, MeanDecreaseGini))
      
    p_imp <- ggplot(imp_df, aes(x = MeanDecreaseGini, y = mz)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = paste0(dataset_name, ": Top 20 Most Important Features"),
           subtitle = "(Predicting segmentation class)",
           x = "Mean Decrease in Gini Impurity", y = "m/z") +
      theme_minimal()
    print(p_imp)
    
    return(rf_model)
  }, error = function(e) {
    warning(sprintf("Classification failed for '%s': %s", dataset_name, e$message))
    return(NULL)
  })
}


#' Construct and plot a co-localization network.
#'
#' @param msi_obj The MSImageSet object.
#' @param dataset_name The name of the dataset.
#' @param threshold_quantile The correlation quantile to use for creating edges.
#' @return An igraph graph object.
construct_network <- function(msi_obj, dataset_name, threshold_quantile = 0.99) {
  if (!loaded_pkgs["igraph"]) return(NULL)
  tryCatch({
    int_mat <- iData(msi_obj)
    corr_mat <- cor(t(int_mat), use = "pairwise.complete.obs")
    diag(corr_mat) <- 0
    
    thr <- quantile(corr_mat[upper.tri(corr_mat)], threshold_quantile, na.rm = TRUE)
    adj_mat <- ifelse(corr_mat >= thr, corr_mat, 0)
    
    g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
    
    igraph::plot(g, vertex.size = 5, vertex.label = NA,
                 edge.width = igraph::E(g)$weight * 5,
                 main = paste0(dataset_name, ": Co-localization Network"))
    return(g)
  }, error = function(e) {
    warning(sprintf("Network construction failed for '%s': %s", dataset_name, e$message))
    return(NULL)
  })
}

## --------------------------------------------------------------------
## Main Execution Loop
## --------------------------------------------------------------------

classifiers <- list()
networks <- list()

purrr::imap(data_list, function(msi_obj, name) {
  message(sprintf("\n--- Starting extra analysis for: %s ---", name))
  
  # Extract features and labels
  feat_res <- extract_features(msi_obj, name, nfeat = 200)
  if (is.null(feat_res)) return(NULL)
  
  # Perform manifold learning (t-SNE/UMAP)
  perform_manifold(feat_res$matrix, feat_res$labels, name)
  
  # Supervised classification and feature importance
  rf_mod <- perform_classification(feat_res$matrix, feat_res$labels, name)
  if (!is.null(rf_mod)) classifiers[[name]] <<- rf_mod
  
  # Co-localization network
  net <- construct_network(msi_obj, name)
  if (!is.null(net)) networks[[name]] <<- net
})

# Save results to the global environment.
assign("msi_classifiers", classifiers, envir = .GlobalEnv)
assign("msi_networks", networks, envir = .GlobalEnv)

message("\nExtra analyses complete.")

## --------------------------------------------------------------------
# INTERPRETATION:
#
# - The t-SNE and UMAP plots help visualize the structure of your data.
#   If the segmentation classes (colors) form distinct islands, it suggests
#   the segments are spectrally unique.
#
# - The Feature Importance plot shows the m/z values the Random Forest
#   model found most useful for telling the segments apart. These are
#   strong candidates for being unique biomarkers for each region.
#
# NEXT STEP:
# Create a new script (e.g., 08_annotate_features.R) to match the
# important m/z values from this script and from the PCA loadings
# against a biological database (e.g., METLIN, HMDB) to hypothesize
# their identities.
## --------------------------------------------------------------------