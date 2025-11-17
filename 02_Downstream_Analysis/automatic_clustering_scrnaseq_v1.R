########## 05/02/2025 ##########

#' @title Automated Clustering Quality Control for scRNA-seq Seurat Objects
#'
#' @description
#' This function performs a fully automated **clustering quality control workflow**
#' for a Seurat object containing multiple clustering resolutions.  
#' It generates UMAPs, barplots, feature plots, boxplots, clustree visualizations,
#' combined QC figures, and a full clustering report, saving all results into a
#' structured sequential directory.  
#' Optionally, the function also carries out cell-cycle–related QC and saves
#' intermediate results.
#'
#' @param SeuratObject A Seurat object that has already undergone clustering at
#' multiple resolutions (i.e., contains metadata columns named
#' \code{"RNA_snn_res.X"} for various resolution values).
#'
#' @param where_to_save Path where all results and plots will be saved.  
#' If \code{NULL} (default), the current working directory is used.
#'
#' @param save_intermediates Logical.  
#' If \code{TRUE} (default), saves the processed Seurat object as an RDS file and
#' reloads it afterward.
#'
#' @param cell_cycle_analysis Logical.  
#' If \code{TRUE} (default), includes cell-cycle scores (S and G2M) in the QC
#' visualizations.
#'
#' @details
#' The workflow includes:
#' \enumerate{
#'   \item Installation and loading of required packages.
#'   \item Loading of custom helper functions.
#'   \item Creation of a dedicated output directory.
#'   \item Iterative analysis across multiple clustering resolutions
#'         (\code{0.1, 0.3, 0.5, 0.7, 1.0}).
#'   \item For each resolution:
#'       \itemize{
#'         \item UMAP visualization of clusters.
#'         \item Barplot showing cluster composition per sample.
#'         \item QC FeaturePlots for indicators such as
#'               \code{nCount_RNA}, \code{nFeature_RNA}, \code{mitoRatio},
#'               and optionally \code{S.Score}, \code{G2M.Score}.
#'         \item Boxplots of metrics across clusters.
#'         \item Combined QC figure summarizing clustering performance.
#'       }
#'   \item Generation of report-level figures:
#'       \itemize{
#'         \item Combined DimPlot across all resolutions.
#'         \item Clustree plot visualizing resolution transitions.
#'         \item UMAP-embedded clustree overlay.
#'         \item Comprehensive report figure combining all QC analyses.
#'       }
#'   \item Optional saving and reloading of processed Seurat object.
#' }
#'
#' @return
#' The input Seurat object, returned without modification unless
#' \code{save_intermediates = TRUE}, in which case the final object is reloaded
#' from an RDS file.  
#' Numerous QC plots and clustering reports are saved to disk in the generated
#' output directories.
#'
#' @seealso
#' \code{\link[Seurat]{DimPlot}},  
#' \code{\link[Seurat]{FeaturePlot}},  
#' \code{\link[clustree]{clustree}},  
#' \code{\link[clustree]{clustree_overlay}}
#'
#' @examples
#' \dontrun{
#' result <- clustering_scrnaseq(
#'     SeuratObject = integrated_seurat,
#'     where_to_save = "results/",
#'     save_intermediates = TRUE,
#'     cell_cycle_analysis = TRUE
#' )
#' }
#'
#' @export

clustering_scrnaseq <- function(SeuratObject, 
                                where_to_save = NULL, 
                                save_intermediates = TRUE,
                                cell_cycle_analysis = TRUE){
  
  # ===================================== #
  # ==== BLOCK 1: Initial Setting Up ====
  # ===================================== #
  # installing packages 
  print('Loading requiered packages')
  list.of.packages = c('ggplot2', 'openxlsx', 'parallel', 
                       'dplyr', 'cowplot', "patchwork", "clustree")
  
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
  
  saving_directory <- create_sequential_dir(path = where_to_save, name = "Clustering")

  # ======================================== #
  # ==== BLOCK 2: Analyzing Resolutions ====
  # ======================================== #
  print_centered_note("Starting Clustering Analysis ")

  print("Starting Clustering Quality Control")
  resolutions = c(0.1, 0.3, 0.5, 0.7, 1)
  resolution_plotlist <- list()
  i <- 1
  for(resolution in resolutions){
    print_centered_note(paste0("Processing Resolution ", resolution))
    resolution_path = create_sequential_dir(saving_directory, 
                                            name = paste0('Clustering_Res_', resolution))
    
    individual_plots <- create_sequential_dir(resolution_path, name = "Individual_Plots")
    
    
    Idents(object = SeuratObject) <- paste0("RNA_snn_res.", resolution)
    ident_index <- grep(paste0("RNA_snn_res.",resolution), colnames(SeuratObject@meta.data))
    
    # One DimPlot for each resolution
    resolution_plot <- DimPlot(SeuratObject,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 6) + NoLegend()+ ggtitle(paste("Resolution ", resolution))

    save_ggplot(resolution_plot, title = paste0("UMAP_Cluster_Res_", resolution),
                folder = individual_plots, width = 2000, height = 3000)


    # ============================================= #
    # ==== BLOCK 3: CLUSTERING QUALITY CONTROL ====
    # ============================================= #
    print_centered_note("Creating QC Clustering Visualizations")
    
    # Barplot of proportion of cell in each cluster by sample
    cat("\n1.Barplot sample cells for each cluster\n")
    
    stacked_barplot <- ggplot(SeuratObject@meta.data) +
      geom_bar(aes(x = SeuratObject@meta.data[[ident_index]], 
                   fill=orig.ident), position=position_fill()) +
      xlab(paste0("Resolution ",resolution)) + 
      theme(legend.position = "bottom", legend.direction = "horizontal")
    
    save_ggplot(stacked_barplot, title = paste0("Barplot_Proportions_Cells_By_Sample_Res_", resolution),
                folder = individual_plots, width = 2000, height = 3000)
    
    # Unwanted sources of variation 
    # Determine metrics to plot present in seurat_integrated@meta.data
    metrics <- c("nCount_RNA", "nFeature_RNA", "mitoRatio")
    if (cell_cycle_analysis) {
      metrics <- c(metrics, "S.Score", "G2M.Score")
    }
    
    print('2.Feature Plot of Sources of Variation')
    
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
    
    # Create FeaturePlots for all metrics
    feature_plots <- lapply(metrics, function(m) generate_feature_plot(SeuratObject, m))
    
    # Create BoxPlots for all metrics
    box_plots <- lapply(metrics, function(m) generate_boxplot(SeuratObject, m, ident_index, resolution))
    
    # Combined plot
    if(cell_cycle_analysis){
      cat("\n3.Combined plot\n")
      combined_plot <- (resolution_plot + stacked_barplot) | (feature_plots[[1]] / feature_plots[[2]] / feature_plots[[3]]) / feature_plots[[4]] / feature_plots[[5]]
      
      combined_plot <- combined_plot + plot_annotation(tag_levels = "A")
      
    }else{
      cat("\n3.Combined plot\n")
      combined_plot <- (resolution_plot + stacked_barplot) | (feature_plots[[1]] / feature_plots[[2]] / feature_plots[[3]]) 
      
      combined_plot <- combined_plot + plot_annotation(tag_levels = "A")
      
    } #else statment key

    save_ggplot(combined_plot, title = paste0("QC_Clustering_Resolution_", resolution), 
                width = 6000, height = 4000, folder = resolution_path)

    resolution_plotlist[[i]] <- resolution_plot
    i <- i+1
    
  } # key for loop resolutions

  # ================================ #
  # ==== BLOCK 4: Report Figure ====
  # ================================ #
  print_centered_note("Creating Clustering Report Plots ")
  individual_plots <- create_sequential_dir(path = saving_directory, 
                                            name = "Individual_Report_Plots")
            
  # Dimplots
  combined_dimplots <- plot_grid(plotlist = resolution_plotlist, nrow = 1)
  
  cat("\nDimPlot of all resolutions\n")
  save_ggplot(plot = combined_dimplots, title = "DimPlot_All_Resolutions", 
              folder = individual_plots,
              width = 6000, height = 3000)
  
  # Clustree plots
  custom_colors <- RColorBrewer::brewer.pal(n = length(resolutions), "Set2")
  
  # umap clustree
  clustree_df <- SeuratObject@meta.data
  clustree_df <- cbind(clustree_df, Embeddings(SeuratObject, reduction = "umap"))
  
  # Crear el Clustree normal con colores personalizados
  clustree <- clustree(clustree_df, prefix = "RNA_snn_res.") + scale_color_manual(values = custom_colors)
  
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
  
  if(save_intermediates){
    print_centered_note("Saving SeuratObject")
    saveRDS(SeuratObject, file = file.path(saving_directory, "Seurat_Clustering_Object.rds"))
    
    rm(SeuratObject, combined_plot, umap_clustree, clustree, clustree_df, resolution_plotlist)
    gc()
    
    SeuratObject <- readRDS(file.path(saving_directory, "Seurat_Clustering_Object.rds"))
  }
  
  print_centered_note("End of Clustering Analysis ")
  
  return(SeuratObject)
} # key of function

