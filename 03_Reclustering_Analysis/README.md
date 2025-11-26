# Single-Cell Analysis Toolkit

### Utilities for Reclustering and Marker Visualization in Seurat

This repository provides a set of R functions designed to simplify
common workflows in single‑cell RNA‑seq analysis using **Seurat**. It
includes:

-   A complete pipeline for **reclustering** selected cell populations
    (`reclustering_function()`).
-   A function for **marker identification and visualization**
    (`plot_markers_function()`).

Both functions automate analysis steps, generate organized outputs, and
produce publication‑quality figures.

------------------------------------------------------------------------

# 📌 Overview

This toolkit provides:

### 🔹 Reclustering Pipeline

`reclustering_function()` enables reclustering of a specific cluster
within an existing Seurat object---useful when a cluster appears
heterogeneous and may contain multiple cellular subtypes.

### 🔹 Marker Visualization Pipeline

`plot_markers_function()` finds the top marker genes for each cluster
and creates multiple visualizations to support cell‑type annotation.

------------------------------------------------------------------------

# 🚀 Features

## 🔥 Reclustering (`reclustering_function()`)

-   Automatic extraction of a selected cluster
-   Full preprocessing: normalization, HVGs, scaling
-   PCA, UMAP, t‑SNE
-   Automatic PC estimation
-   Clustering at multiple resolutions
-   Extensive QC and clustering visualizations
-   Automatic saving of plots and reclustered objects
-   Reproducible directory structure

------------------------------------------------------------------------

## 🎯 Marker Visualization (`plot_markers_function()`)

-   Identifies markers using **FindAllMarkers()**
-   Filters and ranks markers
-   Generates:
    -   DimPlot with labels
    -   DotPlot of top markers
    -   ViolinPlot per cluster
-   Combines and saves visualization panels
-   Exports a sorted marker list to Excel
-   Creates sequential directories to keep outputs organized

------------------------------------------------------------------------

# 📂 Directory Structure

## Reclustering Output

    Reclustering_<cluster_name>/
    ├── Clustering_Res_0.1/
    │   ├── 01_Individual_Plots/
    │   └── 02_QC_Clustering_Resolution_0.1.png
    ├── Clustering_Res_0.3/
    ├── Clustering_Res_0.5/
    ├── Clustering_Res_0.7/
    ├── Clustering_Res_1/
    ├── 01_Reclustering_Dimensionality_Reduction.png
    ├── 06_Clustree_Report.png
    └── Reclustered_<cluster_name>.rds

## Marker Visualization Output

    Plot_Seurat_Markers_001/
    ├── Clustering_Markers_Dataframe.xlsx
    └── Combined_Markers.png

Each run creates a new sequential directory for reproducibility.

------------------------------------------------------------------------

# 🧬 Usage

## Reclustering

``` r
reclustered_obj <- reclustering_function(
    seurat_object = my_seurat_object,
    cluster_to_reanalyze = "3",
    where_to_save = "/path/to/output/"
)
```

### Integrating Results Back Into the Original Object

``` r
new_subclusters <- reclustered_obj$subcluster
seurat_object$subcluster[names(new_subclusters)] <- new_subclusters
```

------------------------------------------------------------------------

## Marker Visualization

``` r
plot_markers_function(
  SeuratObject = my_seurat_object,
  where_to_save = "~/results",
  title = "Seurat_Markers",
  width = 6000,
  height = 4000
)
```

------------------------------------------------------------------------

# 📦 Dependencies

Both functions verify and automatically install required packages:

-   Seurat
-   ggplot2
-   dplyr
-   cowplot
-   patchwork
-   openxlsx
-   clustree (for reclustering)

They also depend on helper scripts:

-   `print_centered_note_v1.R`
-   `create_sequential_dir.R`
-   `Automate_Saving_ggplots.R`
-   `automate_saving_dataframes_xlsx_format.R`

------------------------------------------------------------------------

# 🤝 Contributing

Suggestions and contributions are welcome!\
Feel free to submit pull requests or open issues to extend the toolkit.

------------------------------------------------------------------------

# ✨ Author

**Raúl --- Bioinformatician**
