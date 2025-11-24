# Reclustering Pipeline for Seurat Objects

### A reproducible workflow for subclustering single-cell RNA-seq data

## 📌 Overview

`reclustering_function()` is an automated and fully documented R
workflow designed for **reclustering a specific cluster** within an
existing Seurat object.\
It performs a complete preprocessing, dimensionality reduction,
clustering, and visualization pipeline on a selected subset of cells,
and saves all outputs in a structured directory.

This tool is particularly useful when a cluster appears heterogeneous
and may contain **multiple cellular subtypes**.

## 🚀 Features

✅   Automatic extraction and processing of the selected cluster\
✅   Full Seurat workflow: normalization → HVGs → scaling → PCA → UMAP →
    tSNE
✅   Optimal PC estimation
✅   Clustering at multiple resolutions
✅   Extensive QC and clustering visualizations
✅   Structured output folder system
✅   Automatic saving of all plots and the reclustered Seurat object
✅   Designed for reproducibility

## 📂 Directory Structure

    Reclustering_<cluster_name>/
    ├── Clustering_Res_0.1/
        ├── 01_Individual_Plots/
        └── 02_QC_Clustering_Resolution_0.1.png
    ├── Clustering_Res_0.3/
    ├── Clustering_Res_0.5/
    ├── Clustering_Res_0.7/
    ├── 01_Reclustering_Dimensionality_Reduction.png
    ├── 06_Clustree_Report.png
    └── Reclustered_<cluster_name>.rds

## 🧬 Usage

``` r
reclustered_obj <- reclustering_function(
    seurat_object = my_seurat_object,
    cluster_to_reanalyze = "3",
    where_to_save = "/path/to/output/"
)
```

## 📝 Integration With the Original Seurat Object

``` r
new_subclusters <- reclustered_obj$subcluster
seurat_object$subcluster[names(new_subclusters)] <- new_subclusters
```

## 📦 Dependencies

The function checks and loads: - Seurat\
- ggplot2\
- dplyr\
- cowplot\
- patchwork\
- clustree\
- openxlsx

Plus custom helper scripts.

## 🤝 Contributing

Contributions and suggestions are welcome!


## ✨ Author

**Raúl -- Bioinformatician**
