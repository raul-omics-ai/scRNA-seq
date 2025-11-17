########## 20/01/2025 ##########

#' @title Automated Integration of Seurat Objects Using Harmony
#'
#' @description
#' This function performs a complete automated workflow to integrate multiple
#' Seurat objects using the **Harmony** algorithm.  
#' It includes: installation/loading of required packages, normalization,
#' variable feature selection, PCA, automatic determination of the number of
#' PCs to integrate, generation of plots (elbowplots, pre/post-integration UMAP),
#' execution of Harmony, clustering, and optional saving of intermediate results.  
#' All outputs (plots, RDS files, etc.) are stored in an automatically generated
#' sequential directory.
#'
#' @param list_seurat A list of Seurat objects to be integrated.  
#' Each element must be a valid Seurat object, and the list should be named
#' because the names are used as \code{cell.ids} during merging.
#'
#' @param where_to_save Path where results will be saved.  
#' If \code{NULL} (default), the current working directory is used.
#'
#' @param save_intermediates Logical.  
#' If \code{TRUE} (default), the integrated Seurat object is saved as an RDS file
#' and subsequently reloaded.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Installs and loads necessary packages.
#'   \item Imports custom helper functions required by the pipeline.
#'   \item Creates a directory to store all results.
#'   \item Merges the input Seurat objects.
#'   \item Normalizes data, selects variable features, scales, and runs PCA.
#'   \item Automatically determines the optimal number of PCs using a heuristic
#'         elbowplot approach.
#'   \item Executes UMAP prior to integration.
#'   \item Runs Harmony to integrate datasets by \code{orig.ident}.
#'   \item Runs UMAP, FindNeighbors, and FindClusters on the integrated reduction.
#'   \item Generates and saves pre- and post-integration visualizations.
#'   \item Optionally saves the integrated object and performs memory cleanup.
#' }
#'
#' @return A Seurat object integrated using Harmony, containing the \code{harmony}
#' reduction, UMAP embeddings, and clustering based on the integrated dimensions.  
#' Multiple figures and files are also written into the output directory.
#'
#' @seealso
#' \code{\link[Seurat]{RunHarmony}},  
#' \code{\link[Seurat]{NormalizeData}},  
#' \code{\link[Seurat]{FindVariableFeatures}},  
#' \code{\link[Seurat]{RunUMAP}},  
#' \code{\link[Seurat]{FindClusters}}
#'
#' @examples
#' \dontrun{
#' integrated_seurat <- harmony_automate_integration(
#'     list_seurat = list(sample1 = seurat1, sample2 = seurat2),
#'     where_to_save = "results/",
#'     save_intermediates = TRUE
#' )
#' }
#'
#' @export
harmony_automate_integration <- function(list_seurat, 
                                         where_to_save = NULL, 
                                         save_intermediates = T) {

  # ============================= #
  # ==== BLOCK 1: Setting Up ====
  # ============================= #
  list.of.packages = c('ggplot2','openxlsx', 'parallel', 
                       'harmony', 'dplyr', 'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)

  # custom functions
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")

  # working directory
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  output_dir <- create_sequential_dir(path = where_to_save, name = "Harmony_Integration")

  # ====================================== #
  # ==== BLOCK 2: Harmony Integration ====
  # ====================================== #
  set.seed(123)
  print_centered_note("Starting Integration Studies ")

  print("Merging Seurat Object")
  SeuratObject <- merge(x = list_seurat[[1]],
                        y = list_seurat[2:length(list_seurat)],
                        merge.data = TRUE, add.cell.ids = names(list_seurat))
  
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

  # ========================================= #
  # ==== BLOCK 3: PREINTEGRATION STUDIES ====
  # ========================================= #
  
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
  
  figure <- preintegration_dimplot | postintegration_dimplot
  figure <- figure + plot_annotation(tag_levels = "A")
  
  save_ggplot(plot = figure, title = "Integration_Dimplot", 
              folder = output_dir, height = 2000, width = 3000)
  
  if(save_intermediates){
    print("Saving Integrated Seurat Object")
    saveRDS(SeuratObject, file = file.path(output_dir, "Integrated_Seurat.rds"))
    
    rm(SeuratObject, postintegration_dimplot, preintegration_dimplot,
       dims_to_int)
    gc()
    
    SeuratObject <- readRDS(file = file.path(output_dir, "Integrated_Seurat.rds"))
  }
  
  print_centered_note("End Of The Integration Step")
  return(SeuratObject)
}
