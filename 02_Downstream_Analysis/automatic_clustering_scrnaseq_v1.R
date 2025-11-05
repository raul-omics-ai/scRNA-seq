########## 05/02/2025 ##########

# ========================================= #
# AUTOMATIC CLUSTERING AND QC OF CLUSTERING #
# ========================================= #

clustering_scrnaseq <- function(SeuratObject, 
                                where_to_save = NULL, 
                                save_intermediates = T,
                                cell_cycle_analysis = T){
  # BLOCK 1: LOADING LIBRARIES AND SETTING UP WORKING DIRECTORY ====
  # installing packages 
  print('Loading requiered packages')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'harmony', 'dplyr', 'cowplot', "patchwork", "clustree")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # creating the saving directory
  if(is.null(where_to_save)){
    where_to_save <- getwd()
  }
  
  saving_directory <- file.path(where_to_save, "06_CLUSTERING")
  dir.create(saving_directory, recursive = T, showWarnings = F)
  
  # header
  separator <- paste0(rep("-", 80), collapse = "")
  
  cat("\n")
  print(separator)
  print("                         STARTING CLUSTERING ANALYSIS                           ")
  print(separator)
  
  # BLOCK 2: ANALYZING RESOLUTIONS ====
  print("Starting Clustering Quality Control")
  resolutions = c(0.1, 0.3, 0.5, 0.7, 1)
  resolution_plotlist <- list()
  i <- 1
  for(resolution in resolutions){
    print(paste0('--------------- PROCESSING RESOLUTION  ', resolution, ' ---------------'))
    resolution_path <- file.path(saving_directory, paste0('Clustering_Res_', resolution))
    dir.create(resolution_path, recursive = T)
    
    Idents(object = SeuratObject) <- paste0("RNA_snn_res.", resolution)
    ident_index <- grep(paste0("RNA_snn_res.",resolution), colnames(SeuratObject@meta.data))
    
    # One DimPlot for each resolution
    resolution_plot <- DimPlot(SeuratObject,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 6) + NoLegend()+ ggtitle(paste("Resolution ", resolution))

    png(filename = file.path(resolution_path, paste0("01_UMAP_Cluster_Res_", resolution, ".png")), 
        res = 300, width = 2000, height = 3000)
    print(resolution_plot)
    dev.off()
    
    # BLOCK 3: CLUSTERING QUALITY CONTROL ====
    ## Barplot of proportion of cell in each cluster by sample
    print("1.Barplot sample cells for each cluster")
    
    stacked_barplot <- ggplot(SeuratObject@meta.data) +
      geom_bar(aes(x = SeuratObject@meta.data[[ident_index]], 
                   fill=orig.ident), position=position_fill()) +
      xlab(paste0("Resolution ",resolution)) + 
      theme(legend.position = "bottom", legend.direction = "horizontal")
    
    png(filename = file.path(resolution_path, paste0("02_Barplot_Proportions_Cells_By_Sample_Res_", resolution, ".png")), res = 300, width = 2000, height = 3000)
    print(stacked_barplot)
    dev.off()
    
    ## Segregración por otras fuentes de variación no interesantes
    # Determine metrics to plot present in seurat_integrated@meta.data
    metrics <- c("nCount_RNA", "nFeature_RNA", "mitoRatio")
    if (cell_cycle_analysis) {
      metrics <- c(metrics, "S.Score", "G2M.Score")
    }
    
    print('2.Feature Plot of Sources of Variation')
    
    # Función para generar FeaturePlot
    generate_feature_plot <- function(seurat_obj, feature) {
      FeaturePlot(seurat_obj, 
                  reduction = "umap", 
                  features = feature,
                  pt.size = 0.4, 
                  order = TRUE,
                  min.cutoff = 'q10',
                  label = TRUE)
    }
    
    # Función para generar BoxPlot
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
    feature_plots <- lapply(metrics, function(m) generate_feature_plot(SeuratObject, m))
    
    # Generar BoxPlots para todas las métricas
    box_plots <- lapply(metrics, function(m) generate_boxplot(SeuratObject, m, ident_index, resolution))
    
    # Combined plot
    if(cell_cycle_analysis){
      print("3.Combined plot")
      combined_plot <- (resolution_plot + stacked_barplot) | (feature_plots[[1]] / feature_plots[[2]] / feature_plots[[3]]) / feature_plots[[4]] / feature_plots[[5]]
      png(filename = file.path(resolution_path, paste0("00_Clustering_QC_Resolution_", resolution, ".png")), res = 300, width = 6000, height = 5000)
      print(combined_plot)
      dev.off()
    }else{
      print("3.Combined plot")
      combined_plot <- (resolution_plot + stacked_barplot) | (feature_plots[[1]] / feature_plots[[2]] / feature_plots[[3]]) 
      png(filename = file.path(resolution_path, paste0("00_Clustering_QC_Resolution_", resolution, ".png")), res = 300, width = 6000, height = 5000)
      print(combined_plot)
      dev.off()
    } #else statment key

    resolution_plotlist[[i]] <- resolution_plot
    i <- i+1
    
  } # key for loop resolutions
  
  # Combinar los FeaturePlots en una fila
  print("Saving DimPlots of all dimensions together")
  png(filename = file.path(saving_directory, "01_DimPlot_All_Resolutions.png"), res = 300, width = 6000, height = 3000)
  print(plot_grid(plotlist = resolution_plotlist, nrow = 1))
  dev.off()
  
  # Clustree plots
  print("Creating Clustree plots")
  # Clustree plots
  custom_colors <- RColorBrewer::brewer.pal(n = length(resolutions), "Set2")
  
  # umap clustree
  clustree_df <- SeuratObject@meta.data
  clustree_df <- cbind(clustree_df, Embeddings(SeuratObject, reduction = "umap"))
  
  # Crear el Clustree normal con colores personalizados
  clustree <- clustree(clustree_df, prefix = "RNA_snn_res.") + scale_color_manual(values = custom_colors)
  
  png(filename = file.path(saving_directory, "02_Clustree_Resolutions.png"), res = 300, width = 3000, height = 3000)
  print(clustree)
  dev.off()
  umap_clustree <- clustree_overlay(clustree_df, prefix = "RNA_snn_res.", x_value = "UMAP_1", y_value = "UMAP_2") + 
    scale_fill_manual(values = rev(custom_colors)) + 
    NoLegend()
  
  png(filename = file.path(saving_directory, "03_UMAP_Clustree_Resolutions.png"), res = 300, width = 3000, height = 3000)
  print(umap_clustree)
  dev.off()
  
  # Combined figure
  combined_plot <- (umap_clustree + clustree) / (plot_grid(plotlist = resolution_plotlist, nrow = 1))
  
  png(filename = file.path(saving_directory, "03_Clustree_Report.png"), res = 300, width = 6000, height = 5000)
  print(combined_plot)
  dev.off()
  
  if(save_intermediates){
    saveRDS(SeuratObject, file = file.path(saving_directory, "Seurat_Clustering_Object.rds"))
    
    rm(SeuratObject, combined_plot, umap_clustree, clustree, clustree_df, resolution_plotlist)
    gc()
    
    SeuratObject <- readRDS(file.path(saving_directory, "Seurat_Clustering_Object.rds"))
  }
  
  cat("\n")
  print(separator)
  print("                              END OF CLUSTERING QC                              ")
  print(separator)
  
  return(SeuratObject)
} # key of function

