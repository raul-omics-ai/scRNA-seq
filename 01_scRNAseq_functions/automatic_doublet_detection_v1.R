########## 03/02/2025 ##########

#######################################################
########## AUTOMATE DOUBLET DETECTION SCRIPT ##########
#######################################################

# This function will be an update of the original part of doublet detectio from my automate scRNAseq
# analysis script. I'll be use as a stencil the script created by BioStatsquid, taht uses two different
# tools for doublet-detection (DoubletFinder and scDblFinder) to get a final consensus to take the 
# decission of which cells will be removed.

# In this time, I'll write a function that will be used in the script to simplify code and make things 
# easiest. 

# This function will split seurat object based on orig.ident column from metadata to split by source
# of each sample. If this column doesn't contain this information, add it prior to run this function.

# REFERENCE: https://github.com/biostatsquid/scripts/blob/main/Analysis/scRNAseq/1_Preprocessing/scRNAseq_doublet_detection.R

# NOTE: The doublet detection step will be carry out before any integration step and in each sample
# individually. First split by sample, then identify doublets, merge all datasets and finaly remove
# all doublets.

automatic_doublet_detection <- function(SeuratObject, 
                                        where_to_save = NULL,
                                        save_intermediates = T){
  # BLOCK 1: Library Loading and set up working directory ####
  # INSTALLING PACKAGES 
  options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
  set.seed(42)
  
  print('Loading packages...')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'dplyr', 'cowplot', "scDblFinder", "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # Setting Up Working Directory
  print("Setting Up the Working Directory...")
  if (is.null(where_to_save)) {
    where_to_save <- getwd()
  }
  output_dir <- file.path(where_to_save, paste0("04_Doublet_Detection"))
  dir.create(output_dir, recursive = TRUE)
  
  # Loading DoubletFinder Functions
  separator <- paste0(rep("-", 80), collapse= "")
  cat("\n")
  print(separator)
  print("                          STARTING DOUBLET DETECTION                           ")
  print(separator)
  
  print("Loading custom Functions")
  source("~/Documentos/09_scripts_R/edit_paramsweep_seurat.R")
  source("~/Documentos/09_scripts_R/edited_doubletFinder.R")
  
  doublet_wb <- createWorkbook()
  
  # BLOCK 2: Splitting Seurat Object by biological sample ####
  if(length(unique(SeuratObject$orig.ident)) == 1){
    stop("Remember that orig.ident column must contain the information of the sample source.")
    } else{
    print("These are the samples that will be doublet-detected: ")
    print(unique(as.character(SeuratObject$orig.ident)))
    }
  
  seurat_sample_list <- SplitObject(SeuratObject, split.by = "orig.ident")
  
  doublet_report <- as.data.frame(matrix(nrow = 0, ncol = 7))
  colnames(doublet_report) <- c('Sample', 'Doublet_Detection', 'Singlets', 'Doublets', "Optimal pK", "Multiplet Rate", "Expected Number of Doublets")
  
  for(sample_origin in 1:length(unique(SeuratObject$orig.ident))){
    print(paste0("Detecting Doublets for Sample: ", names(seurat_sample_list)[sample_origin]))
    seu_sample_subset <- seurat_sample_list[[sample_origin]] 
    
    sample_name <- names(seurat_sample_list)[sample_origin]
  
    doublet_report[sample_origin, 1] <- sample_name
    doublet_report[sample_origin, 2] <- 'DoubletFinder'
    # Pre-process seurat object with standard seurat workflow --- 
    sample <- NormalizeData(seu_sample_subset)
    sample <- FindVariableFeatures(sample)
    sample <- ScaleData(sample)
    sample <- RunPCA(sample, nfeatures.print = 10)
    
    # Find significant PCs
    stdv <- sample[["pca"]]@stdev
    percent_stdv <- (stdv/sum(stdv)) * 100
    cumulative <- cumsum(percent_stdv)
    co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
    co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                         percent_stdv[2:length(percent_stdv)]) > 0.1), 
                decreasing = T)[1] + 1
    min_pc <- min(co1, co2)
    # Create a dataframe with values
    plot_df <- data.frame(pct = percent_stdv, 
                          cumu = cumulative, 
                          rank = 1:length(percent_stdv))
    
    # Elbow plot to visualize 
    png(file.path(output_dir, paste0("1_Elbowplot_", sample_name, ".png")), res = 300, width = 2000, height = 2000)
    print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > min_pc)) + 
            geom_text() + 
            geom_vline(xintercept = 90, color = "grey") + 
            geom_hline(yintercept = min(min_pc[min_pc > 5]), color = "grey") +
            theme_bw())
    dev.off()
    # Finish pre-processing with min_pc
    sample <- RunUMAP(sample, dims = 1:min_pc)
    sample <- FindNeighbors(object = sample, dims = 1:min_pc) 
    sample <- FindClusters(object = sample, resolution = 0.1)
    
    # pK identification (no ground-truth) 
    # introduces artificial doublets in varying props, merges with real data set and 
    # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
    # provides a list of the proportion of artificial nearest neighbours for varying
    # combinations of the pN and pK
    sweep_list <- edit_paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
    sweep_stats <- summarizeSweep(sweep_list)
    bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
    # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
    optimal.pk <- bcmvn %>% 
      dplyr::filter(BCmetric == max(BCmetric)) %>%
      dplyr::select(pK)
    optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
    
    doublet_report[sample_origin, 5] <- optimal.pk
    
    # Homotypic doublet proportion estimate
    annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
    homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
    
    # Get the multiplet rate if not provided
    print('Estimating multiplet rate from cells in dataset')
      
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
    
    doublet_report[sample_origin, 6] <- multiplet_rate
    
    nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
    
    doublet_report[sample_origin, 7] <- nExp.poi.adj
    
    # run DoubletFinder
    print('Running DoubletFinder')
    sample <- edited_doubletFinder(seu = sample, 
                            PCs = 1:min_pc, 
                            pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                            nExp = nExp.poi.adj) # number of expected real doublets
    
    # change name of metadata column with Singlet/Doublet information
    colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
    
    print('Creating opt pK plot...')
    png(filename = file.path(output_dir, paste0("2_Optimal_pK_sample_", sample_name, ".png")), res = 300, width = 3000, height = 3000)
    print(ggplot(bcmvn, mapping= aes(pK, BCmetric, group = 1)) +
            geom_point() +
            geom_line())
    dev.off()
    
    # Creando el gráfico combinado para de detección de doublets
    
    # DimPlot para visualizar los doublets en UMAP
    plot1 <- DimPlot(sample, group.by = 'doublet_finder', reduction = 'umap', pt.size = 4) + 
      ggtitle("DimPlot Doublet Detection") + theme(plot.title = element_text(hjust = 0.5)) 
    
    # Violin plots de métricas de calidad
    plot2 <- VlnPlot(sample, group.by = 'doublet_finder', features = "nFeature_RNA", pt.size = 0) + 
      ggtitle("nFeature_RNA") + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
    
    plot3 <- VlnPlot(sample, group.by = 'doublet_finder', features = "nCount_RNA", pt.size = 0) + 
      ggtitle("nCount_RNA") + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
    
    plot4 <- VlnPlot(sample, group.by = 'doublet_finder', features = "mitoRatio", pt.size = 0) + 
      ggtitle("percent_mt") + theme(plot.title = element_text(hjust = 0.5))
    
    # Combinar gráficos en una sola figura
    final_plot <- (plot1 / (plot2 | plot3 | plot4)) + 
      patchwork::plot_annotation(title = paste0("Doublet Detection - ", sample_name))
    
    png(filename = file.path(output_dir, paste0("Doublet_Detection_", sample_name, ".png")), 
        res = 300, width = 4000, height = 4000)
    print(final_plot)
    dev.off()
    
    print(paste0("Doublets detected for sample ", sample_name))
    seurat_sample_list[[sample_origin]] <- sample
    
    doublet_report[sample_origin, 3] <- nrow(sample@meta.data[which(sample@meta.data$doublet_finder == 'Singlet'),])
    doublet_report[sample_origin, 4] <- nrow(sample@meta.data[which(sample@meta.data$doublet_finder == 'Doublet'),])
    
  }# Key for for loop
  
  addWorksheet(doublet_wb, sheetName = "Doublet_Report")
  writeDataTable(doublet_wb, sheet = "Doublet_Report", x = doublet_report)
  
  saveWorkbook(doublet_wb, file = file.path(output_dir, "Doublet_Report.xlsx"), overwrite = T)
  
  print("Merging all dataset into one Seurat Object")
  SeuratObject <- merge(x = seurat_sample_list[[1]],
                                 y = seurat_sample_list[2:length(seurat_sample_list)],
                                 merge.data = TRUE, add.cell.ids = names(seurat_sample_list))
  
  print("Removing Doublets from the dataset")
  SeuratObject <- subset(SeuratObject,
                     subset = doublet_finder == 'Singlet')
  
  if(save_intermediates){
    print("Saving Merged Seurat Object")
    saveRDS(object = SeuratObject, file.path(output_dir, "Seurat_Doublet_Free_Object.rds"))
    
    rm(SeuratObject, doublet_wb, doublet_report, multiplet_rates_10x, seurat_sample_list)
    gc()
    
    SeuratObject <- readRDS(file = file.path(output_dir, "Seurat_Doublet_Free_Object.rds"))
  }
  cat("\n")
  print(separator)
  print("                           END OF DOUBLET DETECTION                             ")
  print(separator)
  
  return(SeuratObject)
} # Key to close function
