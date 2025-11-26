########## 24/11/2025 ##########

# ======================================== #
# ==== Marker identification function ====
# ======================================== #

# Esta función lo que pretende es crear visualizaciones acerca de los marcadores de cada uno de los clusters
# para poder tomar la decisión de anotar el tipo celular que se corresponda

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




