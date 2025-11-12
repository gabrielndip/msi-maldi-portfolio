## 08_annotate_features.R
##
## Biological Annotation of Significant m/z Features
##
## This script takes the most significant m/z values identified in the
## PCA (loadings) and Random Forest (feature importance) analyses and
## attempts to assign putative molecular identities to them by querying
## the METLIN public database.
##
## NOTE: This script requires an internet connection to function.
## Annotation is a complex, research-grade problem. The results from this
## script should be considered **hypotheses** that require further
## validation, not definitive identifications.

suppressPackageStartupMessages({
  library(tidyverse)
  library(httr)
  library(jsonlite)
})

## --------------------------------------------------------------------
## Configuration
## --------------------------------------------------------------------

# The number of top features to extract from PCA and RF analyses.
n_top_features <- 10

# The mass tolerance for database searching in parts-per-million (ppm).
# This is a critical parameter. A value between 5-10 is reasonable for
# modern TOF instruments. Adjust based on your instrument's mass accuracy.
mass_tolerance_ppm <- 10

# The adducts to search for. This example focuses on common positive
# ion mode adducts. Add or change these based on your experimental setup.
# Examples: "[M+H]+", "[M+Na]+", "[M+K]+"
adducts_to_search <- c("[M+H]+", "[M+Na]+")

## --------------------------------------------------------------------
## 1. Gather Significant m/z Features
## --------------------------------------------------------------------

message("Gathering significant m/z features from previous analyses...")

# --- From PCA Loadings (Script 06) ---
pca_features <- list()
if (exists("msi_pca") && is.list(msi_pca) && length(msi_pca) > 0) {
  pca_features <- purrr::map(msi_pca, function(pca_res) {
    loadings_df <- as.data.frame(loadings(pca_res))
    loadings_df$mz <- Cardinal::mz(pca_res)
    
    # Get top features for each PC
    purrr::map(1:ncol(loadings(pca_res)), function(i) {
      pc_col <- paste0("PC", i)
      loadings_df %>%
        arrange(desc(abs(.data[[pc_col]]))) %>%
        slice_head(n = n_top_features) %>%
        pull(mz)
    }) %>% unlist() %>% unique()
  })
}

# --- From Random Forest Importance (Script 07) ---
rf_features <- list()
if (exists("msi_classifiers") && is.list(msi_classifiers) && length(msi_classifiers) > 0) {
  rf_features <- purrr::map(msi_classifiers, function(rf_mod) {
    if (!is.null(rf_mod)) {
      randomForest::importance(rf_mod) %>%
        as.data.frame() %>%
        rownames_to_column("mz") %>%
        arrange(desc(MeanDecreaseGini)) %>%
        slice_head(n = n_top_features) %>%
        pull(mz) %>%
        as.numeric()
    }
  })
}

# Combine and get unique m/z values across all analyses
significant_mz <- c(unlist(pca_features), unlist(rf_features)) %>%
  unique() %>%
  sort()

if (length(significant_mz) == 0) {
  stop("No significant m/z features found. Please run 06_pca_analysis.R and 07_extra_analysis.R first.")
}

message(sprintf("Found %d unique significant m/z values to annotate.", length(significant_mz)))

## --------------------------------------------------------------------
## 2. Annotate Features using METLIN API
## --------------------------------------------------------------------

#' Query METLIN for a single m/z value.
#'
#' @param mz The m/z value to query.
#' @param tolerance_ppm The search tolerance in ppm.
#' @param adduct The adduct to search (e.g., "[M+H]+").
#' @return A tibble of putative annotations, or NULL if none found.
query_metlin <- function(mz, tolerance_ppm, adduct) {
  # The METLIN API endpoint is no longer publicly documented for this type of query.
  # This function provides a placeholder structure for how one would perform such a query.
  # For a real workflow, you would use a dedicated R package like 'MetaboAnnotation'
  # or a web service with a documented API.
  
  # This is a simulated function.
  # In a real scenario, you would make an HTTP request here.
  # Example of what the request might look like:
  # url <- paste0("https://metlin.scripps.edu/rest/v1/name/m/z/", mz, "/tolerance/", tolerance_ppm, "/adduct/", adduct)
  # response <- httr::GET(url)
  # content <- httr::content(response, "text", encoding = "UTF-8")
  # results <- jsonlite::fromJSON(content)
  
  # Since we cannot make live requests, we return an empty tibble
  # to demonstrate the structure.
  return(tibble(
    query_mz = numeric(),
    adduct = character(),
    metlin_id = character(),
    name = character(),
    formula = character(),
    exact_mass = numeric(),
    ppm_error = numeric()
  ))
}

message("Starting annotation... (This is a demonstration and will not perform live web requests)")

# This section is for demonstration purposes.
# The `query_metlin` function is a placeholder.
# To perform real annotation, you would replace it with a call to a real annotation package or API.

all_annotations <- purrr::map_dfr(significant_mz, function(mz) {
  purrr::map_dfr(adducts_to_search, function(adduct) {
    # In a real run, this function would make a web request.
    # query_metlin(mz, mass_tolerance_ppm, adduct)
    
    # For this demo, we just print a message.
    cat(sprintf("  (Demo) Querying for m/z: %.4f with adduct: %s\n", mz, adduct))
    
    # Return an empty tibble to show the desired structure
    tibble()
  })
})


## --------------------------------------------------------------------
## 3. Display Results
## --------------------------------------------------------------------

if (nrow(all_annotations) > 0) {
  message("Annotation complete. Found potential matches:")
  print(all_annotations)
} else {
  message("\nAnnotation demonstration complete. No live queries were performed.")
  message("To perform real annotation, you would need to replace the placeholder `query_metlin` function.")
  message("Consider using a Bioconductor package like 'MetaboAnnotation' for robust, real-world annotation workflows.")
}

# Save the list of significant m/z values for manual searching.
significant_features_df <- tibble(mz = significant_mz)
write_csv(significant_features_df, "significant_mz_features.csv")

message("\nA file 'significant_mz_features.csv' has been created with the list of m/z values for manual database searching.")

## --------------------------------------------------------------------
# INTERPRETATION:
#
# The output table would show potential molecular formulas and names for your
# significant m/z values.
#
# - **ppm_error**: This is the difference between your measured m/z and the
#   theoretical m/z of a database entry. Lower is better.
# - **Multiple Hits**: It is common for a single m/z to have multiple
#   potential matches (isomers).
#
# These results are starting points for further investigation. You can use
# this information to search literature or design follow-up experiments
# (e.g., MS/MS fragmentation) to confirm identities.
## --------------------------------------------------------------------
