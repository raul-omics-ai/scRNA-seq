########## 04/02/2025 ##########

# ========================================== #
# AUTOMATE COMPLETE SCRNASEQ ANALYSIS SCRIPT #
# ========================================== #

# This scripts will host the new automate scRNAseq function that will be used from scratch, improving the
# code, make the function more efficient and easy to read. 

# The main difference with the previous script is in this new function, only will be passed one scRNA-seq
# that contains the experiment I want to analyze, not a list of scRNA-seq, so it can be used with lapply
# to analyze multiple scRNA-seq experiments at onece. 

# I'm not going to introduce the cell cycle analysis at begining, but I hope that in later upgrates of the
# function will be added. Also this function will create new and better visualizations, and a better performance
# in saving results, figures and so.

# The same as above with the integration step that will be performed, at firts, out of the function.

automatic_scrnaseq_analysis <- function(SeuratObject, where_to_save = NULL,
                                        cell_cycle_analysis = TRUE, 
                                        title = "scRNAseq_Analysis", 
                                        specie = "hsa",
                                        doublet_detection = TRUE,
                                        integration = TRUE, 
                                        save_intermediates_files = TRUE,
                                        filter_nUMIs = 500,
                                        filter_nGenes = 300,
                                        filter_mitoRatio = 0.2,
                                        filter_log10GenesPerUMI = 0){
  # BLOCK 1: LOADING LIBRARIES, CUSTOM FUNCTIONS AND CREATING NEW DIRECTORIES ====
  # installing packages 
  print('Instalando los paquetes necesarios para el análisis...')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'harmony', 'dplyr', 'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # loading custom functions
  source("~/Documentos/09_scripts_R/automatic_QC_scrnaseq.R")
  source("~/Documentos/09_scripts_R/automate_standard_seurat_processing_v1.R")
  source("~/Documentos/09_scripts_R/automatic_doublet_detection_v1.R")
  source("~/Documentos/09_scripts_R/harmony_automate_integration_v1.R")
  source("~/Documentos/09_scripts_R/automatic_clustering_scrnaseq_v1.R")
  # creating the main directory
  if(is.null(where_to_save)){
    where_to_save <- getwd()}
  
  main_directory <- file.path(where_to_save, paste0(title, "_", Sys.Date()))
  dir.create(main_directory, recursive = T, showWarnings = F)
  
  # BLOCK 2: INITIAL QC ====
  SeuratObject <- automatic_qc_scrnaseq(SeuratObject = SeuratObject, where_to_save = main_directory, specie = specie)
  
  # BLOCK 3: FILTERING ====
  separator <- paste0(rep("-", 80), collapse = "")
  cat("\n")
  print(separator)
  print("                   STARTING WITH FILTERING LOW QUALITY CELLS                    ")
  print(separator)
  
  dir.create(path = file.path(main_directory, "02_FILTERING"), recursive = T, showWarnings = F)
  # Crear un excel con los parámetros de filtro
  wb <- createWorkbook()
  print('Creating report of the filtering')
  
  addWorksheet(wb, sheetName = "Filtering Parameters")
  filtering_parameters <- data.frame(matrix(nrow=4, ncol = 2))
  colnames(filtering_parameters) <- c('QC_metric', 'Filter')
  filtering_parameters$QC_metric <- c('nUMIs (Valor mínimo)', 'nGenes (Valor mínimo)', 'mitoRatio (Valor máximo)', 'log10GenesPerUMI (Valor mínimo)')
  filtering_parameters$Filter <- c(filter_nUMIs, filter_nGenes, filter_mitoRatio, filter_log10GenesPerUMI)
  
  writeData(wb, sheet = "Filtering Parameters", x = filtering_parameters)
  # filtering SeuratObject
  print('Filtering Seurat Object')
  SeuratObject <- subset(SeuratObject, subset = nCount_RNA >= filter_nUMIs &
             nFeature_RNA >= filter_nGenes &
             mitoRatio <= filter_mitoRatio &
             log10GenesPerUMI >= filter_log10GenesPerUMI)
  
  # Filtering report
  filtering_report <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(filtering_report) <- c('Sample', 'QC_Metric', 'Prefiltering', 'PostFiltering')
  
  # Especificar las métricas de calidad
  qc_metrics <- c('N cells', 'Median N Genes', 'Median N UMIs', 'Median mitoRatio', 'Median log10Genes/UMI')
  for(sample in unique(SeuratObject$orig.ident)){
    temp_report <- data.frame(
      Sample = rep(sample, length(qc_metrics)),
      QC_Metric = qc_metrics,
      Prefiltering = NA,
      PostFiltering = NA)
    
    # Ncounts
    temp_report[1, 'Prefiltering'] <- length(SeuratObject@meta.data$orig.ident == sample)
    temp_report[1, 'PostFiltering'] <- length(SeuratObject@meta.data$orig.ident == sample)
    
    # N Genes
    temp_report[3, 'Prefiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "nFeature_RNA"])
    temp_report[3, 'PostFiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "nFeature_RNA"])
    
    # N UMIs
    temp_report[2, 'Prefiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "nCount_RNA"])
    temp_report[2, 'PostFiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "nCount_RNA"])
    
    # mitoRatio
    temp_report[4, 'Prefiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "mitoRatio"])
    temp_report[4, 'PostFiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "mitoRatio"])
    
    # log10Genes/UMI (Complexity)
    temp_report[5, 'Prefiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "log10GenesPerUMI"])
    temp_report[5, 'PostFiltering'] <- median(SeuratObject@meta.data[which(SeuratObject@meta.data$orig.ident == sample), "log10GenesPerUMI"])
    
    # Agregar el dataframe temporal al dataframe principal
    filtering_report <- rbind(filtering_report, temp_report)
  } # Key for loop
  
  filtering_report[,c(3,4)] <- apply(filtering_report[, c(3,4)], 2, round, digits = 2)  
  addWorksheet(wb, sheet = "Filtering Report")
  writeDataTable(wb, sheet = "Filtering Report", x = filtering_report)
  saveWorkbook(wb, file = file.path(main_directory, "02_FILTERING/Filtering_Report.xlsx"))
  
  if(save_intermediates_files){
    print('Saving Seurat Object Post Filtering')
    saveRDS(SeuratObject, file.path(main_directory, "02_FILTERING", "Seurat_Post_filtering.rds"))
    
    rm(SeuratObject, filtering_report,temp_report, sample, filter_nUMIs, filter_nGenes, filter_mitoRatio, filter_log10GenesPerUMI, filtering_parameters, wb)
    gc()
    
    SeuratObject <- readRDS(file.path(main_directory, "02_FILTERING", "Seurat_Post_filtering.rds"))
  }
  
  SeuratObject <- automatic_qc_scrnaseq(SeuratObject = SeuratObject, where_to_save = file.path(main_directory, "02_FILTERING"))
  
  
  
  # BLOCK 3: STANDARD SEURAT PIPELINE ====
  SeuratObject <- automatic_standard_seurat(SeuratObject = SeuratObject, 
                                            where_to_save = file.path(main_directory), 
                                            save_intermediate_files = save_intermediates_files, 
                                            cell_cycle_analysis = cell_cycle_analysis)
  
  # BLOCK 4: DOUBLET DETECTION ====
  if(doublet_detection){
    SeuratObject <- automatic_doublet_detection(SeuratObject = SeuratObject, 
                                                where_to_save = file.path(main_directory), 
                                                save_intermediates = save_intermediates_files)
  }
  
  
  # BLOCK 5: INTEGRATION STUTDY ====
  if(integration){
    SeuratObject <- harmony_automate_integration(SeuratObject = SeuratObject, 
                                                 where_to_save = file.path(main_directory), 
                                                 save_intermediates = save_intermediates_files)
  }
  
  
  # BLOCK 6: CLUSTERING AND CLUSTERING QC
  SeuratObject <- clustering_scrnaseq(SeuratObject = SeuratObject, 
                                      where_to_save = file.path(main_directory), 
                                      save_intermediates = save_intermediates_files, 
                                      cell_cycle_analysis =cell_cycle_analysis)
  return(SeuratObject)
  
} # Key for the end of the function
