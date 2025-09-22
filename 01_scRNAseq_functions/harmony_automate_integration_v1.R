########## 20/01/2025 ##########

## Script automatizado para el análisis de integración con harmony de un objeto Seurat:

harmony_automate_integration <- function(SeuratObject, where_to_save = NULL, save_intermediates = T) {
  # BLOCK 1: PREPARING WORKING DIRECTORY ==== 
  print('Instalando los paquetes necesarios para el análisis...')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'harmony', 'dplyr', 'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  print("Setting Up the Working Directory")
  if (is.null(where_to_save)) {
    where_to_save <- getwd()
  }
  output_dir <- file.path(where_to_save, paste0("05_INTEGRATION"))
  dir.create(output_dir, recursive = TRUE)
  
  # header
  separator <- paste0(rep("-", 80), collapse = "")
  cat("\n")
  print(separator)
  print("                          STARTING INTEGRATION STUDIES                          ")
  print(separator) 
  
  # BLOCK 2: PROCESSING MERGED SEURAT OBJECT ====
  print("Processing merged Seurat Object")
  SeuratObject <- SeuratObject %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% 
    ScaleData() %>%
    RunPCA()
  
  print("Creating Elbowplots")
  pct <- SeuratObject[["pca"]]@stdev / sum(SeuratObject[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  pcs <- ifelse(!is.na(co1) & !is.na(co2), min(co1, co2), max(length(pct), na.rm = TRUE))
  
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  elbow_plot_path <- file.path(output_dir, "Merged_dataset_Elbowplot.png")
  png(elbow_plot_path, res = 300, width = 2000, height = 2000)
  print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
          geom_text() + 
          geom_vline(xintercept = 90, color = "grey") + 
          geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
          theme_bw())
  dev.off()
  
  dims_to_int <- pcs
  SeuratObject <- RunUMAP(SeuratObject, dims = 1:dims_to_int)
  
  # BLOCK 3: PREINTEGRATION STUDIES ====
  preintegration_dimplot <- DimPlot(SeuratObject, reduction = "umap", group.by = "orig.ident") + 
    ggtitle("PreIntegration") +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  
  print(paste0("A total of  ", dims_to_int, " dimensions will be used to integrate all datasets"))
  
  print("Starting the integration step with Harmony")
  SeuratObject <- RunHarmony(SeuratObject, 
                                  group.by.vars = "orig.ident", 
                                  reduction = "pca", 
                                  reduction.save = "harmony")
  
  SeuratObject <- SeuratObject %>%
    RunUMAP(reduction = "harmony", dims = 1:dims_to_int) %>%
    FindNeighbors(reduction = "harmony", dims = 1:dims_to_int) %>%
    FindClusters(resolution = c(0.1, 0.3, 0.5, 0.7, 1.0))
  
  print("Integration finished")
  print("Creating the visualizations of integration")
  
  postintegration_dimplot <- DimPlot(SeuratObject, group.by = "orig.ident", reduction = "umap") + 
          ggtitle("PostIntegration") + NoLegend()
  
  integrated_umap_path <- file.path(output_dir, "01_Integration_Dimplots.png")
  
  png(integrated_umap_path, res = 300, width = 3000, height = 2000)
  print(preintegration_dimplot | postintegration_dimplot)
  dev.off()
  
  if(save_intermediates){
    print("Saving Integrated Seurat Object")
    saveRDS(SeuratObject, file = file.path(output_dir, "Seurat_Integrado.rds"))
    
    rm(SeuratObject, postintegration_dimplot, integrated_umap_path, preintegration_dimplot,
       dims_to_int, elbow_plot_path)
    gc()
    
    SeuratObject <- readRDS(file = file.path(output_dir, "Seurat_Integrado.rds"))
  }
  cat("\n")
  print(separator)
  print("                            END OF INTEGRATION STEP                             ")
  print(separator)
  return(SeuratObject)
}
