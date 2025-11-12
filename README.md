# MSI-MALDI Portfolio: Brain Tissue Analysis

## 1. Project Overview

This project provides a complete workflow for the analysis of Matrix-Assisted Laser Desorption/Ionization (MALDI) Mass Spectrometry Imaging (MSI) data, demonstrated on a brain tissue dataset. The primary goal is to process raw MSI data, identify distinct molecular regions within the tissue, and extract key spectral features differentiating these regions.

The analysis is conducted in R, primarily using the `Cardinal` package for MSI data handling and analysis. The workflow is structured as a series of sequential scripts, each performing a distinct step from data loading to advanced analysis.

## 2. Table of Contents
- [1. Project Overview](#1-project-overview)
- [2. Table of Contents](#2-table-of-contents)
- [3. Data Description](#3-data-description)
- [4. Workflow](#4-workflow)
- [5. Getting Started](#5-getting-started)
  - [5.1. Prerequisites](#51-prerequisites)
  - [5.2. Installation](#52-installation)
  - [5.3. Running the Analysis](#53-running-the-analysis)
- [6. Script-by-Script Breakdown](#6-script-by-script-breakdown)
- [7. Results](#7-results)

## 3. Data Description

The `data/` directory contains the MSI datasets. The data is in the standard imzML format, which consists of two files per dataset:
- **`.imzML`**: An XML file containing metadata about the experiment (e.g., scan settings, spatial coordinates).
- **`.ibd`**: A binary file containing the mass spectral intensity data.

The datasets included are:
- `Brain01_Bregma-1-46_centroid`: The primary brain tissue sample.
- `cal-nonormalization`, `control-nonormalization`, `png1-nonormalization`, `png2-nonormalization`: Additional datasets, likely for calibration or as controls.

## 4. Workflow

The analysis pipeline is executed through a series of R scripts located in the `scripts/` directory. The general workflow is as follows:

1.  **Load Data**: Read the `.imzML` and `.ibd` files into R.
2.  **Preprocess**: Normalize intensities, pick, align, and filter spectral peaks.
3.  **Explore**: Generate quality control plots (TIC images, spectra) to inspect the data.
4.  **Crop to ROI**: Automatically identify and crop the image to the tissue area (Region of Interest).
5.  **Segmentation**: Group pixels into spatially-aware clusters using Spatial Shrunken Centroids (SSC).
6.  **Dimensionality Reduction**: Perform PCA to identify major sources of variation.
7.  **Advanced Analysis**: Further explore the data with t-SNE, UMAP, and co-localization networks.

## 5. Getting Started

### 5.1. Prerequisites

- R (version 4.0 or later recommended)
- RStudio (recommended)

### 5.2. Installation

Open R and run the following command to install the required packages from CRAN and Bioconductor:

```R
# Install core packages from CRAN
install.packages(c("tidyverse", "Rtsne", "uwot", "randomForest", "igraph"))

# Install Bioconductor manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Cardinal from Bioconductor
BiocManager::install("Cardinal")
```

### 5.3. Running the Analysis

The analysis can be run in two main ways:

1.  **Interactive Execution (Recommended)**:
    Open the `msi_maldi_portfolio_workflow.qmd` file in RStudio. This Quarto document provides a narrative walkthrough of the analysis. You can run each code chunk sequentially to see the outputs and visualizations directly.

2.  **Running Scripts Sequentially**:
    Alternatively, you can run the R scripts in the `scripts/` directory in numerical order:
    - `01_load_data.R`
    - `02_preprocess_data.R`
    - `03_explore_data.R`
    - `04_roi_crop.R`
    - `05_ssc_segmentation.R`
    - `06_pca_analysis.R`
    - `07_extra_analysis.R`

## 6. Script-by-Script Breakdown

- **`01_load_data.R`**: Loads all `.imzML` datasets from the `data/` directory.
- **`02_preprocess_data.R`**: Applies TIC normalization, peak picking, alignment, and filtering.
- **`03_explore_data.R`**: Generates diagnostic plots (TIC images, ion images, spectra) for quality control.
- **`04_roi_crop.R`**: Automatically crops the datasets to the tissue region to remove empty space.
- **`05_ssc_segmentation.R`**: Performs spatial segmentation to identify anatomically distinct regions based on their mass spectra.
- **`06_pca_analysis.R`**: Runs Principal Component Analysis (PCA) to find the main patterns of spectral variation across the tissue.
- **`07_extra_analysis.R`**: Demonstrates advanced techniques like t-SNE/UMAP for visualization, and co-localization network analysis.

## 7. Results

The primary output of this project is the rendered HTML report:
- **`msi_maldi_portfolio_workflow.html`**: This file contains the full analysis, including code, explanatory text, and all generated figures (e.g., PCA score plots, segmentation maps). Open this file in a web browser to view the complete results.

Figures generated during the analysis are saved in the `msi_maldi_portfolio_workflow_files/figure-html/` directory.
