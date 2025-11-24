########## 24/11/2025 ##########

# =============================================================== #
# ==== Reclustering Pipeline for a Selected Seurat Cluster  ==== #
# ============================================================== #

#' Reclustering Pipeline for a Selected Seurat Cluster
#'
#' @description
#' This function performs a full reclustering workflow on a selected cluster
#' from a Seurat object. It subsets the specified cluster, re-runs the essential
#' Seurat preprocessing steps (normalization, scaling, PCA, UMAP, tSNE),
#' performs clustering across multiple resolutions, generates extensive QC and
#' clustering visualizations, saves all output plots and intermediate files,
#' and finally returns the reclustered Seurat object.
#'
#' The function is designed as a self-contained reclustering module for
#' identifying substructure within a given cluster (e.g., potential subtypes).
#'
#' @details
#' The workflow consists of four major blocks:
#'
#' \enumerate{
#'   \item{\strong{Initial Setup:}}{
#'     Loads required packages, checks user input, sources helper functions,
#'     and creates the output directory using a sequential naming system.
#'   }
#'
#'   \item{\strong{Subset & Full Seurat Processing:}}{
#'     Subsets the specified cluster and runs:
#'     \code{NormalizeData()}, \code{FindVariableFeatures()},
#'     \code{ScaleData()}, \code{RunPCA()}, optimal PC estimation,
#'     \code{RunUMAP()}, and \code{RunTSNE()}.
#'     Dimensionality reduction plots are generated and saved.
#'   }
#'
#'   \item{\strong{Clustering Analysis Across Resolutions:}}{
#'     Runs \code{FindNeighbors()} and \code{FindClusters()} for a predefined
#'     set of resolutions. For each resolution, the function generates:
#'     \itemize{
#'       \item UMAP cluster plots
#'       \item Cell proportion barplots per sample
#'       \item FeaturePlots for QC metrics (nCount_RNA, nFeature_RNA, mitoRatio)
#'       \item Boxplots for QC metrics
#'       \item Combined QC/clustering panels
#'     }
#'     All plots are automatically saved.
#'   }
#'
#'   \item{\strong{Clustering Report:}}{
#'     Generates clustree visualizations and UMAP-overlaid clustree plots.
#'     A combined report is also created and saved.
#'   }
#' }
#'
#' The final reclustered Seurat object is saved as an RDS file in the designated
#' directory.
#'
#' @param seurat_object A Seurat object. Must contain valid identities
#'   (\code{Idents}) that include \code{cluster_to_reanalyze}.
#'
#' @param cluster_to_reanalyze Character. The identity (cluster name) 
#'   to subset and recluster.
#'
#' @param where_to_save Character or \code{NULL}. Path to the directory where
#'   output files will be saved. If \code{NULL}, the current working directory is used.
#'
#' @return
#' A reclustered Seurat object with all dimensionality reduction embeddings,
#' cluster assignments, and metadata corresponding to the reclustering workflow.
#'
#' @section Side Effects:
#' This function produces several side effects:
#' \itemize{
#'   \item Creates a folder titled \code{Reclustering_<clustername>} using
#'         a sequential directory generator.
#'   \item Saves numerous ggplot2 figures (UMAP, tSNE, PCA, QC, clustree).
#'   \item Saves the final reclustered Seurat object as an \code{.rds} file.
#' }
#'
#' @note
#' This function is intended for reclustering of a subset of cells and does not
#' merge the reclustered object back into the original object. To reintegrate
#' subcluster annotations, extract the metadata from the returned object and
#' append it to the original Seurat object.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[Seurat]{FindClusters}}
#'   \item \code{\link[Seurat]{RunUMAP}}
#'   \item \code{\link[clustree]{clustree}}
#' }
#'
#' @examples
#' \dontrun{
#' # Reclustering cluster "3" and saving inside the working directory
#' reclust_obj <- reclustering_function(
#'     seurat_object = my_seurat,
#'     cluster_to_reanalyze = "3"
#' )
#'
#' # Saving to a specific output folder
#' reclust_obj <- reclustering_function(
#'     seurat_object = my_seurat,
#'     cluster_to_reanalyze = "Macrophages",
#'     where_to_save = "/home/user/results/"
#' )
#' }
#'
#' @export

reclustering_function <- function(seurat_object,
                                  cluster_to_reanalyze,
                                  where_to_save = NULL){
  # ==================================================================== #
  # ==== BLOCK 1: Set up functions, libraries and working directory ====
  # ==================================================================== #
  # checkpoint
  if(!cluster_to_reanalyze %in% Idents(seurat_object)){
    stop("Please review the cluster name to reanalyze. \nThere may be a spelling error.")
  }
  
  # loading packages
  list.of.packages = c('ggplot2', 'openxlsx', "Seurat",'dplyr', 
                       'cowplot', "patchwork", "clustree")
  
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
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  saving_directory <- create_sequential_dir(path = where_to_save, 
                                            name = paste0("Reclustering_", cluster_to_reanalyze))
  
  # ============================================================================== #
  # ==== BLOCK 2: Complete Seurat Processing Pipeline in Subset Seurat Object ====
  # ============================================================================== #
  print_centered_note("Subsetting the main Seurat Object")
  
  reclustering_object <- subset(x = seurat_object, idents = cluster_to_reanalyze)
  
  print_centered_note("Continue with the Seurat Pipeline")
  reclustering_object <- reclustering_object %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
  
  print("Computing the optimal number of PCs to select")
  # Determine percent of variation associated with each PC
  pct <- reclustering_object[["pca"]]@stdev / sum(reclustering_object[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  print("Performing Non-linear Dimensionality Reduction Step")
  dims_per_sample <- pcs
  print(paste0('A total of ', pcs, ' PCs were selected for UMAP and TSNE'))
  reclustering_object <- RunUMAP(reclustering_object, dims = 1:dims_per_sample)
  print('RunUMAP OK')
  reclustering_object <- RunTSNE(reclustering_object, dims = 1:dims_per_sample)
  print('RunTSNE OK')
  
  print_centered_note("Creating Visualizations for Dimensonality Reduction")
  plot_list <- list()
  for(reduction in c("pca", "umap", "tsne")){
    plot <-  DimPlot(reclustering_object, reduction = reduction) + 
      ggtitle(toupper(reduction)) +
      theme(legend.position = "bottom")
    
    plot_list[[reduction]] <- plot
  }
  plot <- plot_grid(plotlist = plot_list, ncol = 3)
  
  save_ggplot(plot, title = "Reclustering_Dimensionality_Reduction", 
              folder = saving_directory,
              dpi = 300, width = 3000, height = 2000)
  
  # ====================================== #
  # ==== BLOCK 3: Clustering Analysis ====
  # ====================================== #
  print_centered_note("Performing Clustering Analysis ")
  
  resolutions = c(0.1, 0.3, 0.5, 0.7, 1)
  
  reclustering_object <- reclustering_object %>%
    FindNeighbors(reduction = "umap", dims = 1:2) %>%
    FindClusters(resolution = resolutions)
  
  resolution_plotlist <- list()
  i <- 1
  for(resolution in resolutions){
    print_centered_note(paste0("Processing Resolution ", resolution))
    
    resolution_path = create_sequential_dir(saving_directory, 
                                            name = paste0('Clustering_Res_', resolution))
    
    individual_plots <- create_sequential_dir(resolution_path, name = "Individual_Plots")
    
    Idents(object = reclustering_object) <- paste0("RNA_snn_res.", resolution)
    ident_index <- grep(paste0("RNA_snn_res.",resolution), colnames(reclustering_object@meta.data))
    
    # One DimPlot for each resolution
    resolution_plot <- DimPlot(reclustering_object,
                               reduction = "umap",
                               label = TRUE,
                               label.size = 6) + NoLegend()+ ggtitle(paste("Resolution ", resolution))
    
    save_ggplot(resolution_plot, title = paste0("UMAP_Cluster_Res_", resolution),
                folder = individual_plots, width = 2000, height = 3000)
  
    print_centered_note("Creating QC Clustering Visualizations")
    
    # Barplot of proportion of cell in each cluster by sample
    cat("\n1.Barplot sample cells for each cluster\n")
    
    stacked_barplot <- ggplot(reclustering_object@meta.data) +
      geom_bar(aes(x = reclustering_object@meta.data[[ident_index]], 
                   fill=orig.ident), position=position_fill()) +
      xlab(paste0("Resolution ",resolution)) + 
      theme(legend.position = "bottom", legend.direction = "horizontal")
    
    save_ggplot(stacked_barplot, title = paste0("Barplot_Proportions_Cells_By_Sample_Res_", resolution),
                folder = individual_plots, width = 2000, height = 3000)
    
    
    # Unwanted sources of variation 
    # Determine metrics to plot present in seurat_integrated@meta.data
    metrics <- c("nCount_RNA", "nFeature_RNA", "mitoRatio")
    
    cat('\n2.Feature Plot of Sources of Variation\n')
    
    # Function to generate FeaturePlots
    generate_feature_plot <- function(seurat_obj, feature) {
      FeaturePlot(seurat_obj, 
                  reduction = "umap", 
                  features = feature,
                  pt.size = 0.4, 
                  order = TRUE,
                  min.cutoff = 'q10',
                  label = TRUE)
    }
    
    # Function to generate BoxPlots
    generate_boxplot <- function(seurat_obj, feature, ident_index, resolution) {
      ggplot(seurat_obj@meta.data) +
        geom_boxplot(aes(x = seurat_obj@meta.data[[ident_index]], 
                         y = seurat_obj@meta.data[[feature]], 
                         fill = factor(seurat_obj@meta.data[[ident_index]]))) +
        NoLegend() +
        xlab(paste0("Resolution ", resolution)) +
        ylab(feature)
    }
    
    # Generar FeaturePlots para todas las métricas
    feature_plots <- lapply(metrics, function(m) generate_feature_plot(reclustering_object, m))
    
    # Generar BoxPlots para todas las métricas
    box_plots <- lapply(metrics, function(m) generate_boxplot(reclustering_object, m, ident_index, resolution))
    
    # Combined plo
    cat("\n3.Combined plot\n")
    combined_plot <- (resolution_plot + stacked_barplot) | (feature_plots[[1]] / feature_plots[[2]] / feature_plots[[3]]) 
    
    combined_plot <- combined_plot + plot_annotation(tag_levels = "A")
    
    save_ggplot(combined_plot, title = paste0("QC_Clustering_Resolution_", resolution), 
                width = 6000, height = 4000, folder = resolution_path)
    
    resolution_plotlist[[i]] <- resolution_plot
    i <- i+1
    } # key for loop resolutions
  
  # =================================================== #
  # ==== BLOCK 4: Clustering Report Visualizations ====
  # =================================================== #
  print_centered_note("Creating the Visualizations for Clustering ")
  
  # Clustree plots
  custom_colors <- RColorBrewer::brewer.pal(n = length(resolutions), "Set2")
  
  # umap clustree
  clustree_df <- reclustering_object@meta.data
  clustree_df <- cbind(clustree_df, Embeddings(reclustering_object, reduction = "umap"))
  
  # Crear el Clustree normal con colores personalizados
  clustree <- clustree(clustree_df, prefix = "RNA_snn_res.") + 
    scale_color_manual(values = custom_colors)
  
  cat("\nClustree Plot\n")
  save_ggplot(plot = clustree, 
              title = "Clustree_Resolutions", 
              folder = individual_plots, 
              width = 3000, height = 3000)
  
  cat("\nUMAP Clustree Plot\n")
  umap_clustree <- clustree_overlay(clustree_df, prefix = "RNA_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2") + 
    scale_fill_manual(values = rev(custom_colors)) + 
    NoLegend()
  
  save_ggplot(plot = umap_clustree, 
              title = "UMAP_Clustree_Resolutions", 
              folder = individual_plots, 
              width = 3000, height = 3000)
  
  # Combined figure
  combined_plot <- (umap_clustree + clustree) / (plot_grid(plotlist = resolution_plotlist, nrow = 1))
  
  combined_plot <- combined_plot + plot_annotation(tag_levels = "A")
  
  cat("\nCombined Figure\n")
  save_ggplot(plot = combined_plot, title = "Clustree_Report", 
              folder = saving_directory,
              width = 6000, height = 5000)
  
  print_centered_note("Saving Reclustered Seurat Object ")
  saveRDS(reclustering_object, 
          file = file.path(saving_directory, paste0("Reclustered_", cluster_to_reanalyze, ".rds")))
  cat("\nReclustered Object Saved\n")
  
  print_centered_note("End of the Reclustering Analysis ")
  return(reclustering_object)
} # main function key
