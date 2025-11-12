<!-- 
------------------------------------------------------------------------
## YYYY-MM-DD: [Brief Title of Today's Work]
------------------------------------------------------------------------

**Goal:** 
- [A short sentence on the main objective for the day. E.g., "Test the effect of different 's' parameters in SSC segmentation."]

**Method / Changes:**
- [List the specific changes made to code or parameters.]
- `scripts/05_ssc_segmentation.R`: Changed `ssc_s` from 15 to `c(5, 10, 20)`.
- [Note any other significant actions, e.g., "Downloaded a new dataset from Zenodo."]

**Observations:**
- [What was the outcome? What did you notice?]
- Segmentation with `s=5` was too noisy.
- `s=20` produced the most visually distinct and anatomically relevant regions.
- The script failed when...

**Next Steps:**
- [What will you do next based on these observations?]
- Proceed with `s=20` for the final segmentation.
- Investigate the co-localization network of the new segments.

-->

---
# Previous Notes
---

## 1. Choose and download a larger dataset

-   A good candidate is the publicly available mouse kidney dataset prepared by Melanie Föll and colleagues and hosted on Zenodo (derived from PRIDE PXD009808). It contains MALDI images of formalin‑fixed murine kidney tissue sections with N‑linked glycans. Three 6‑µm sections were imaged; two were treated with PNGase F to release N‑glycans, while a third contained control and calibrant regions. The data were acquired at 100µm spatial resolution and processed to create concise imzML files (m/z range 1250‑2310, resampled at 0.1 Da, with TIC normalization) zenodo.org

-   The Zenodo record includes the combined dataset (all_files.imzml/.ibd, \~16.8 MB + 466 MB) as well as individual “control” and “treated1” sections

-   ***Download from PRIDE repo:***

    ***\$ wget -r -nH -np --cut-dirs=7 <ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2018/11/PXD009808/>***

## 2  Load and preprocess the larger data

1.  Set the file path: Download the .imzML and corresponding .ibd files for the chosen dataset. In your R script, change the path in readMSIData() to the new file:

*Example: mouse kidney dataset (combined sections)*

*data_path \<- "path/to/all_files.imzML"*

*mse \<- readMSIData(data_path, mass.range = c(1250, 2310), resolution = 0.01, units = "mz")*

2.  Limit the mass range (optional): If the dataset covers a broad m/z range, load only the range of interest using the mass.range argument or downsample using resolution. This dramatically reduces memory usage.
3.  Crop to a region of interest (ROI): Before expensive processing, subset the pixel grid to a relevant tissue region. For the kidney dataset, you might remove the calibrant area by selecting x/y ranges:

*roi_pixels \<- pixels(mse, x \>= 500 & x \<= 2000, y \>= 300 & y \<= 1600)*

*mse_roi \<- mse[features(mse), roi_pixels]*

4.  Adjust preprocessing parameters:
    -   Normalization: retain TIC normalization (method="tic") or experiment with median normalization if TIC varies strongly.
    -   Peak picking: For high‑resolution data, specify units="ppm" and a narrower tolerance; e.g. peakPick(method="mad", SNR=6) followed by peakAlign(tolerance=5, units="ppm").
    -   Filtering: Increase freq.min (e.g. 0.02 or 0.05) to remove rare peaks and speed up downstream analyses.
5.  Process: Apply process() after the preprocessing chain to commit changes. Save intermediates (saveRDS()) so you can resume later without re‑running the entire pipeline.

## 3  Extend the analysis beyond the original tutorial

Once the larger dataset is loaded and preprocessed, you can perform deeper or alternative analyses that were not in the original tutorial:

### 3.1 Dimensionality reduction and visualization

-   **UMAP / t‑SNE**: Use `umap` or `Rtsne` on the matrix of picked peaks to visualize nonlinear patterns in the data. For example:

    ```         
    library(Rtsne)
    tsne_res <- Rtsne(t(intensityMatrix(mse_proc)), perplexity=30)
    plot(tsne_res$Y, col = label_vector, pch=16)
    ```

-   **Non‑negative matrix factorization (NMF)**: `NMF` can extract additive components corresponding to co-localized ion profiles.

-   **Independent component analysis (ICA)**: Identify statistically independent sources of variation.

### 3.2 Spatially aware segmentation

-   **Alternative clustering**: Besides Spatial‑Shrunken‑Centroids (SSC), try:

    -   **K‑means or hierarchical clustering** on the principal component scores or on selected peaks.

    -   **PhenoGraph** or density‑based clustering (DBSCAN) for irregular shapes.

    -   **Graph cut or watershed** algorithms on ion images to detect anatomical structures.

-   **SPUTNIK**: The R package *SPUTNIK* filters peaks based on spatial autocorrelation and co‑occurrence. Use it to remove isotopes/noise and then re‑cluster.

-   **Morphological segmentation**: Build a tissue mask from the total‑ion image (e.g. thresholding, `EBImage`) and propagate labels into cluster results to focus on biological tissue.

### 3.3 Supervised and comparative analyses

-   **Supervised classification**: If you have labels (e.g. treated vs control kidney sections), train a Random Forest, Support Vector Machine, or PLS‑DA model on the spectral features to identify discriminative ions. Cardinal’s `sda()` function implements shrinkage discriminant analysis.

-   **Cross‑sample integration**: Combine multiple tissue sections (e.g., the control and treated kidney files[[zenodo.org]{.underline}](https://zenodo.org/records/2628280#:~:text=We%20reduced%20the%20m%2Fz%20range,All%20processing%20steps%20were%20performed){alt="https://zenodo.org/records/2628280#:~:text=We%20reduced%20the%20m%2Fz%20range,All%20processing%20steps%20were%20performed"}) using `combine()` from Cardinal or manual concatenation. After batch normalization, perform PCA or clustering to see how samples differ.

-   **Differential ion analysis**: Use `limma` or generalized linear models to test each m/z feature for differential abundance between regions or conditions, adjusting for multiple testing.

### 3.4 Biological annotation and network analysis

-   **Metabolite/Glycan annotation**: Use the provided LC‑MS/MS table (Glycan_IDs) in the kidney dataset[[zenodo.org]{.underline}](https://zenodo.org/records/2628280#:~:text=We%20processed%20the%20original%20imzML,All%20processing%20steps%20were%20performed){alt="https://zenodo.org/records/2628280#:~:text=We%20processed%20the%20original%20imzML,All%20processing%20steps%20were%20performed"} to assign identities to discriminative m/z features. You can integrate with databases (e.g. GlyTouCan) using `MSnbase` or custom scripts.

-   **Co‑localization networks**: Compute pairwise Pearson/Spearman correlations between ion images; build a network where nodes are m/z features and edges represent strong spatial correlation. Community detection can reveal co‑localized metabolic modules.

-   **Pathway mapping**: After annotation, map glycans or lipids to pathways (e.g. glycosylation, glycerophospholipid metabolism) and interpret spatial differences.

### 3.5 Visualization enhancements

-   **Composite ion images**: Overlay three or more ions using RGB channels to highlight co‑localization patterns; adjust contrast and scaling consistently.

-   **Interactive exploration**: Use `plotly` to make interactive spectra or spatial heatmaps that allow zooming into m/z ranges or pixel subsets.

-   **Automated reports**: Wrap your adapted workflow in an `.Rmd` or `quarto` document that generates an HTML report with key QC plots, segmentation maps and lists of top features.

------------------------------------------------------------------------

## 4  Documenting your personalised workflow

Throughout the adaptation, comment your code extensively to explain why parameters were chosen, what each step accomplishes, and how the new analyses address biological questions. Store scripts and reports in a version‑controlled repository so you can show interviewers or collaborators your thought process and proficiency in MSI analysis.

By choosing a larger, real‑world dataset (such as the murine kidney glycan dataset, [[zenodo.org]{.underline}](https://zenodo.org/records/2628280#:~:text=Three%206%C2%B5m%20sections%20of%20formalin,TOF%2FTOF%20instrument){alt="https://zenodo.org/records/2628280#:~:text=Three%206%C2%B5m%20sections%20of%20formalin,TOF%2FTOF%20instrument"}[,]{.underline} or the urinary bladder phospholipid dataset, [[proteomecentral.proteomexchange.org]{.underline}](https://proteomecentral.proteomexchange.org/cgi/GetDataset#:~:text=Title%20Mass%20spectrometry%20imaging%20of,19){alt="https://proteomecentral.proteomexchange.org/cgi/GetDataset#:~:text=Title%20Mass%20spectrometry%20imaging%20of,19"}), and implementing advanced segmentation, dimension reduction, supervised learning and biological annotation, you will move beyond the original tutorial and demonstrate deeper expertise in MALDI‑MSI data science.