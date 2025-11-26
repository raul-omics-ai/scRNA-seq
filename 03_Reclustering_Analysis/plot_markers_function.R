########## 24/11/2025 ##########

#' Plot Marker Genes for Seurat Clusters
#'
#' This function identifies marker genes for each cluster within a Seurat object
#' and generates several visualizations (DimPlot, DotPlot, ViolinPlot), combining them into
#' a single figure. It also exports a sorted marker table and automatically saves all
#' generated outputs in a sequential directory.
#'
#' @description
#' The function computes positive markers for each cluster using \code{FindAllMarkers},
#' extracts the top markers, and generates publication-ready plots to assist in deciding the
#' appropriate cell type annotation for each cluster.
#'
#' Several helper functions are sourced internally:
#'  - \code{print_centered_note()}
#'  - \code{create_sequential_dir()}
#'  - \code{save_ggplot()}
#'  - \code{save_dataframe()}
#'
#' These must be available in the user's environment.
#'
#' @param SeuratObject A valid Seurat object containing clustering results.
#' @param where_to_save Character string specifying the directory in which results will be saved.
#'   If \code{NULL}, the current working directory is used. Defaults to \code{NULL}.
#' @param title A character string defining the name of the main output directory.
#'   Defaults to \code{"Plot_Seurat_Markers"}.
#' @param width Numeric. Width of the combined marker plot when saved. Defaults to \code{6000}.
#' @param height Numeric. Height of the combined marker plot when saved. Defaults to \code{4000}.
#'
#' @return
#' The function invisibly returns the path to the directory where all outputs were saved.
#' In addition, it writes:
#'  - An Excel file with sorted cluster markers.
#'  - A combined figure with DimPlot, DotPlot, and ViolinPlot.
#'
#' @details
#' The function automatically installs missing packages, loads required dependencies,
#' finds all positive markers, extracts the top two genes per cluster, creates visualizations,
#' and saves all outputs in a structured directory.
#'
#' @examples
#' \dontrun{
#' plot_markers_function(
#'   SeuratObject = my_seurat,
#'   where_to_save = "~/results",
#'   title = "My_Marker_Plots"
#' )
#' }
#'
#' @export

plot_markers_function <- function(SeuratObject, 
                                  where_to_save = NULL,
                                  title = "Plot_Seurat_Markers",
                                  width = 6000, height = 4000){
  
  # ======================================= #
  # ==== BLOCK 1: Setting up function ====  #
  # ======================================= #
  # loading packages
  list.of.packages = c('ggplot2', 'openxlsx', "Seurat",'dplyr', 
                       'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # custom functions
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  
  # working directory
  print_centered_note("Creating the Saving Results Directory")
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  saving_directory <- create_sequential_dir(path = where_to_save, 
                                            name = title)
  # =================================== #
  # ==== BLOCK 2: Find all markers ====
  # =================================== #
  print_centered_note("Finding All Markers for Clusters ")
  markers <- FindAllMarkers(SeuratObject, 
                            only.pos = TRUE)
  
  markers_sorted <- markers %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 1) %>%
    #slice_max(order_by = avg_log2FC, with_ties = FALSE) %>%
    ungroup() %>%
    as.data.frame()
  
  cat("\nSaving Marker List\n")
  save_dataframe(markers_sorted, title = "Clustering_Markers_Dataframe", folder = saving_directory)
  
  # select 2 markers for each cluster
  top2_genes <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
    slice_head(n = 2) %>%
    pull(gene) %>%
    unique()
  
  # ================================= #
  # ==== BLOCK 3: Visualizations ====
  # ================================= #
  print_centered_note("Creating the Figure of the Markers ")
  cat("\nDimPlot\n")
  dimplot <- DimPlot(SeuratObject, label = T, repel = T,
                     pt.size = 8,   label.size = 8) + NoLegend()
  
  cat("\nDotPlot\n")
  dotplot <- DotPlot(object = SeuratObject, features = top2_genes)
  
  cat("\nViolinPlot\n")
  violinplot <- VlnPlot(SeuratObject, features = top2_genes)
  
  cat("\nCombining Plots\n")
  combined <- (dotplot / violinplot) | dimplot
  combined <- combined + plot_annotation(tag_levels = "A")
  
  cat("\nSaving Combining Plot\n")
  save_ggplot(combined, title = "Combined_Markers", 
              folder = saving_directory, width = width, height = height)
  
  print_centered_note("End of the function")
} # main function key




