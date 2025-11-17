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

automatic_scrnaseq_analysis <- function(list_scr, 
                                        where_to_save = NULL,
                                        cell_cycle_analysis = TRUE, 
                                        title = "scRNAseq_Analysis", 
                                        specie = "hsa",
                                        doublet_detection = TRUE,
                                        integration = TRUE, 
                                        save_intermediates_files = TRUE){
  
  # =================================================================================== #
  # ==== BLOCK 1: LOADING LIBRARIES, CUSTOM FUNCTIONS AND CREATING NEW DIRECTORIES ====
  # =================================================================================== #   
  # installing package
  print('Setting Up Directory and Loading Packages')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'harmony', 'dplyr', 'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # loading custom functions
  source("~/Documentos/09_scripts_R/automatic_QC_scrnaseq.R")
  source("~/Documentos/09_scripts_R/automate_standard_seurat_processing.R")
  source("~/Documentos/09_scripts_R/automatic_doublet_detection.R")
  source("~/Documentos/09_scripts_R/harmony_automate_integration.R")
  source("~/Documentos/09_scripts_R/automatic_clustering_scrnaseq.R")
  
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  
  # creating the main directory
  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  main_directory <- create_sequential_dir(path = where_to_save, 
                                          name = paste0(title, "_", Sys.Date()))

  
  # ========================================== #
  # ==== BLOCK 2: INITIAL Quality Control ====
  # ========================================== #
  list_scr <- automatic_qc_scrnaseq(list_srn = list_scr, 
                                    where_to_save = main_directory, 
                                    specie = specie)
  
  # ============================================= #
  # === BLOCK 3: LOW QUALITY CELLS FILTERING ====
  # ============================================= #
  print_centered_note("Starting With Filtering Step To Remove Low Quality Cells ")

  filtering_dir <- create_sequential_dir(path = main_directory, name = "Filtering")

  # Creating a filtering parameter report
  wb <- createWorkbook()
  print('Creating report of the filtering')

  filter_scrnaseq_function<- function(One_seurat, wb = wb) {
    sample <- unique(One_seurat$orig.ident)
    cat(paste0("\nFiltering Setting for Sample ", sample, "\n"))
    cat("\nPress ENTER for default parameters\n")
    
    # Interactive thresholds
    filter_nUMIs <- as.numeric(readline(prompt = "Minimum number of UMIs per cell (default 500): "))
    filter_nGenes <- as.numeric(readline(prompt = "Minimum number of genes per cell (default 300): "))
    filter_mitoRatio <- as.numeric(readline(prompt = "Maximum mitoRatio (default 0.2): "))
    filter_log10GenesPerUMI <- as.numeric(readline(prompt = "Minimum complexity (log10GenesPerUMI, default 0.8): "))
    
    # If user just presses Enter, use defaults
    if (is.na(filter_nUMIs)) filter_nUMIs <- 500
    if (is.na(filter_nGenes)) filter_nGenes <- 300
    if (is.na(filter_mitoRatio)) filter_mitoRatio <- 0.2
    if (is.na(filter_log10GenesPerUMI)) filter_log10GenesPerUMI <- 0.8
    
    # Create a table with filtering parameters
    filtering_parameters <- data.frame(
      QC_metric = c('nUMIs (Min Value)', 'nGenes (Min Value)', 
                    'mitoRatio (Max Value)', 'log10GenesPerUMI (Min Value)'),
      Filter = c(filter_nUMIs, filter_nGenes, filter_mitoRatio, filter_log10GenesPerUMI)
    )
    
    # Save to Excel if workbook is provided
    if (!is.null(wb)) {
      addWorksheet(wb, sheetName = paste0(sample, "_params"))
      writeData(wb, sheet = paste0(sample, "_params"), x = filtering_parameters)
    }
    
    seurat_filtered <- subset(
      One_seurat,
      subset = nCount_RNA >= filter_nUMIs &
        nFeature_RNA >= filter_nGenes &
        mitoRatio <= filter_mitoRatio &
        log10GenesPerUMI >= filter_log10GenesPerUMI
    )

    return(seurat_filtered)
  }

  # Filtering Seurat
  filtered_scr_list <- lapply(list_scr, filter_scrnaseq_function, wb)

  # Filtering report
  qc_metrics <- c("Number of Cells", "Median nUMI", "Median nGenes", "Median mitoRatio", "Median log10GenesPerUMI")
  filtering_report <- data.frame()

  # Iterate through all samples in the list
  for (sample_name in names(list_scr)) {
    cat("Processing sample:", sample_name, "\n")
    
    Sort_SeuratObject <- list_scr[[sample_name]]        # Before filtering
    Seurat_filt <- filtered_scr_list[[sample_name]]    # After filtering
    
    # Create temporary dataframe for the sample
    temp_report <- data.frame(
      Sample = rep(sample_name, length(qc_metrics)),
      QC_Metric = qc_metrics,
      Prefiltering = NA,
      PostFiltering = NA
    )
    
    # Fill in the metrics
    # 1. Number of Cells
    temp_report[1, 'Prefiltering'] <- ncol(Sort_SeuratObject)
    temp_report[1, 'PostFiltering'] <- ncol(Seurat_filt)
    
    # 2. Median nUMIs
    temp_report[2, 'Prefiltering'] <- median(Sort_SeuratObject@meta.data$nCount_RNA)
    temp_report[2, 'PostFiltering'] <- median(Seurat_filt@meta.data$nCount_RNA)
    
    # 3. Median nGenes
    temp_report[3, 'Prefiltering'] <- median(Sort_SeuratObject@meta.data$nFeature_RNA)
    temp_report[3, 'PostFiltering'] <- median(Seurat_filt@meta.data$nFeature_RNA)
    
    # 4. Median mitoRatio
    temp_report[4, 'Prefiltering'] <- median(Sort_SeuratObject@meta.data$mitoRatio)
    temp_report[4, 'PostFiltering'] <- median(Seurat_filt@meta.data$mitoRatio)
    
    # 5. Median log10GenesPerUMI
    temp_report[5, 'Prefiltering'] <- median(Sort_SeuratObject@meta.data$log10GenesPerUMI)
    temp_report[5, 'PostFiltering'] <- median(Seurat_filt@meta.data$log10GenesPerUMI)
    
    # Add to main dataframe
    filtering_report <- rbind(filtering_report, temp_report)
  }
  
  # Round numeric values
  filtering_report[, c("Prefiltering", "PostFiltering")] <- 
    apply(filtering_report[, c("Prefiltering", "PostFiltering")], 2, round, digits = 2)
  
  addWorksheet(wb, sheet = "Filtering Report")
  writeDataTable(wb, sheet = "Filtering Report", x = filtering_report)
  saveWorkbook(wb, file = file.path(filtering_dir, "Filtering_Report.xlsx"))

  list_scr <- filtered_scr_list
  
  if(save_intermediates_files){
    print('Saving Seurat Object Post Filtering')
    saveRDS(list_scr, file.path(filtering_dir, "Seurat_Post_filtering.rds"))
    
    rm(Sort_SeuratObject, Seurat_filt, filtering_report,temp_report, sample_name, wb)
    gc()
    
    list_scr <- readRDS(file.path(filtering_dir, "Seurat_Post_filtering.rds"))
  }

  # QC of the filtered data
  list_scr <- automatic_qc_scrnaseq(list_srn = list_scr, 
                                    where_to_save = filtering_dir,
                                   specie = specie)

  # =========================================== #
  # ==== BLOCK 3: STANDARD SEURAT PIPELINE ====
  # =========================================== #
  list_scr <- automatic_standard_seurat(list_seurat = list_scr, 
                                            where_to_save = main_directory, 
                                            save_intermediate_files = save_intermediates_files,
                                            specie = specie,
                                            cell_cycle_analysis = cell_cycle_analysis)
  # ==================================== #
  # ==== BLOCK 4: DOUBLET DETECTION ====
  # ==================================== #
  if(doublet_detection){
    list_scr <- automatic_doublet_detection(list_seurat = list_scr, 
                                                where_to_save = main_directory, 
                                                save_intermediates = save_intermediates_files)
  }
  
  # ===================================== #
  # ==== BLOCK 5: INTEGRATION STUTDY ====
  # ===================================== #
  if(integration){
    SeuratObject <- harmony_automate_integration(list_seurat = list_scr, 
                                                 where_to_save = main_directory, 
                                                 save_intermediates = save_intermediates_files)
  }
  
  # =============================================== #
  # ==== BLOCK 6: CLUSTERING AND CLUSTERING QC ====
  # =============================================== #
  SeuratObject <- clustering_scrnaseq(SeuratObject = SeuratObject, 
                                      where_to_save = main_directory, 
                                      save_intermediates = save_intermediates_files, 
                                      cell_cycle_analysis = cell_cycle_analysis)
  end <- Sys.time()
  execution_time <- end - start
  print_centered_note(paste0("Execution time: ", execution_time, " s"))
  return(SeuratObject)
} # Key for the end of the function
