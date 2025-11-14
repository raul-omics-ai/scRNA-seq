########## 03/02/2025 ##########

#' Automatic Doublet Detection Pipeline for Multiple Seurat Objects
#'
#' Executes a fully automated pipeline for doublet detection in single-cell
#' RNA-seq data using **DoubletFinder** and additional helper utilities.
#'
#' The pipeline includes:
#' \itemize{
#'   \item Standard Seurat preprocessing (normalization, variable feature
#'         selection, scaling, PCA).
#'   \item Automatic detection of the optimal number of PCs.
#'   \item Automated pK parameter search using \code{edit_paramSweep()} and
#'         \code{find.pK()}.
#'   \item Automatic estimation of the multiplet rate based on 10X Genomics
#'         recovery tables.
#'   \item Computation of the adjusted number of expected doublets.
#'   \item Execution of DoubletFinder via \code{edited_doubletFinder()}.
#'   \item Automatic generation and saving of diagnostic plots and reports.
#' }
#'
#' The function processes multiple Seurat objects, generates reports and plots
#' per sample, removes identified doublets, and returns an updated list of
#' doublet-free Seurat objects.
#'
#' @param list_seurat A list of \link[Seurat]{Seurat} objects to process. Each
#'        element must contain the metadata column \code{orig.ident}.
#' @param where_to_save Path where output folders will be created.  
#'        If \code{NULL} (default), the current working directory is used.
#' @param save_intermediates Logical. If \code{TRUE} (default), intermediate
#'        objects such as the doublet-free Seurat list will be saved.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Automatically creates a sequential directory named
#'         \code{Doublet_Detection}.
#'   \item For each Seurat object:
#'     \itemize{
#'       \item Adjusts cell names to avoid duplicates.
#'       \item Runs a standard Seurat preprocessing workflow.
#'       \item Determines the optimal number of PCs using empirical criteria.
#'       \item Estimates pK and multiplet rate.
#'       \item Runs DoubletFinder.
#'       \item Generates diagnostic plots: elbow plot, pK plot, DimPlot and
#'             violin plots.
#'       \item Saves all outputs in structured subfolders.
#'       \item Filters out detected doublets and updates the Seurat object.
#'     }
#'   \item Produces an XLSX summary report with detection statistics per sample.
#' }
#'
#' The following helper scripts must be available:
#' \code{print_centered_note_v1.R},  
#' \code{Automate_Saving_ggplots.R},  
#' \code{automate_saving_dataframes_xlsx_format.R},  
#' \code{create_sequential_dir.R},  
#' \code{edit_paramsweep_seurat.R},  
#' \code{edited_doubletFinder.R}.
#'
#' @return Returns a list of Seurat objects where detected doublets have been
#'         removed.  
#'         If \code{save_intermediates = TRUE}, results are also saved to the
#'         output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load a list of Seurat objects
#' seurat_list <- list(sample1, sample2, sample3)
#'
#' # Run the automated doublet detection workflow
#' seurat_clean <- automatic_doublet_detection(
#'     list_seurat = seurat_list,
#'     where_to_save = "results/doublets/"
#' )
#' }
#' REFERENCE: https://github.com/biostatsquid/scripts/blob/main/Analysis/scRNAseq/1_Preprocessing/scRNAseq_doublet_detection.R

automatic_doublet_detection <- function(list_seurat, 
                                        where_to_save = NULL,
                                        save_intermediates = T){
  # =============================================================== #
  # ==== BLOCK 1: Library Loading and set up working directory ====
  # =============================================================== #
  # installing packages 
  options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
  set.seed(42)
  
  print('Loading packages...')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'dplyr', 'cowplot', "scDblFinder", "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)

  # loading custom functions
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/edit_paramsweep_seurat.R")
  source("~/Documentos/09_scripts_R/edited_doubletFinder.R")
  
  # setting up working directory
  print_centered_note("Setting Up the Working Directory ")  if (is.null(where_to_save)) {

  where_to_save <- ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  output_dir <- create_sequential_dir(path = where_to_save, name = "Doublet_Detection")

  # ======================================== #  
  # ==== BLOCK 2: Running DoubletFinder ==== #
  # ======================================== #  
  print_centered_note("Starting Doublet Detection ")  
  
  doublet_report <- as.data.frame(matrix(nrow = length(list_seurat), ncol = 7))
  colnames(doublet_report) <- c('Sample', 'Doublet_Detection', 'Singlets', 'Doublets',
                                "Optimal pK", "Multiplet Rate", "Expected Number of Doublets")
  parametrization_folder <- create_sequential_dir(output_dir, name = "Parametrization_DoubletFinder")
  
  for (i in seq_along(list_seurat)) {
    sample_name <- unique(list_seurat[[i]]$orig.ident)
    
    new_names <- paste0(sample_name, "_", colnames(list_seurat[[i]]))
    list_seurat[[i]] <- RenameCells(list_seurat[[i]], new.names = new_names)
  }  
  for(i in seq_along(list_seurat)){
    seu_sample <- list_seurat[[i]]
    sample <- unique(seu_sample$orig.ident)
    print_centered_note(paste0("Detecting Doublets for Sample: ", sample))
    
    doublet_report[i, 1] <- as.character(sample)
    doublet_report[i, 2] <- 'DoubletFinder'
    
    # Pre-process seurat object with standard seurat workflow --- 
    seu_sample <- seu_sample %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA(nfeatures.print = 10)
      
    # Find significant PCs
    stdv <- seu_sample[["pca"]]@stdev
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
    elbowplot <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > min_pc)) + 
      geom_text() + 
      geom_vline(xintercept = 90, color = "grey") + 
      geom_hline(yintercept = min(min_pc[min_pc > 5]), color = "grey") +
      theme_bw()
    
    save_ggplot(plot = elbowplot, title = paste0("Elbowplot_", sample), folder = parametrization_folder,
                width = 2000, height = 2000)
    
    # Finish pre-processing with min_pc
    seu_sample <- seu_sample %>%
      RunUMAP(dims = 1:min_pc) %>%
      FindNeighbors(dims = 1:min_pc) %>%
      FindClusters(resolution = 0.1)
    
    # Diagnostics: check cell names / metadata alignment
    cat("\n")
    cat("== Sample:", sample, "==\n")
    cat("ncol counts (if present):", if(!is.null(seu_sample@assays$RNA@counts)) ncol(seu_sample@assays$RNA@counts) else "no counts", "\n")
    cat("length(colnames):", length(colnames(seu_sample)), "\n")
    cat("length(meta rows):", nrow(seu_sample@meta.data), "\n")
    print(head(colnames(seu_sample), 10))
    print(head(rownames(seu_sample@meta.data), 10))
    cat("Any duplicated colnames?:", any(duplicated(colnames(seu_sample))), "\n")
    cat("Any NA colnames?:", any(is.na(colnames(seu_sample))), "\n")
    cat("All colnames == rownames(meta)?", all(colnames(seu_sample) == rownames(seu_sample@meta.data)), "\n")
    cat("\n")

    # pK identification (no ground-truth) 
    # introduces artificial doublets in varying props, merges with real data set and 
    # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
    # provides a list of the proportion of artificial nearest neighbours for varying
    # combinations of the pN and pK
    sweep_list <- edit_paramSweep(seu_sample, PCs = 1:min_pc, sct = FALSE)   
    sweep_stats <- summarizeSweep(sweep_list)
    bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
    
    # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
    optimal.pk <- bcmvn %>% 
      dplyr::filter(BCmetric == max(BCmetric)) %>%
      dplyr::select(pK)
    optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
    
    doublet_report[i, 5] <- optimal.pk
    
    # Homotypic doublet proportion estimate
    annotations <- seu_sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
    homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
    
    # Get the multiplet rate if not provided
    print_centered_note('Estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    
    multiplet_rate <- multiplet_rates_10x %>% 
      dplyr::filter(Recovered_cells < nrow(seu_sample@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
    
    doublet_report[i, 6] <- multiplet_rate
    
    nExp.poi <- round(multiplet_rate * nrow(seu_sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
    
    doublet_report[i, 7] <- nExp.poi.adj
    
    # run DoubletFinder
    print_centered_note('Running DoubletFinder')
    seu_sample <- edited_doubletFinder(seu = seu_sample, 
                                       PCs = 1:min_pc, 
                                       pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                                       nExp = nExp.poi.adj) # number of expected real doublets
    
    # change name of metadata column with Singlet/Doublet information
    colnames(seu_sample@meta.data)[grepl('DF.classifications.*', colnames(seu_sample@meta.data))] <- "doublet_finder"
    
    print('Creating opt pK plot...')
    
    pk_plot <- ggplot(bcmvn, mapping= aes(pK, BCmetric, group = 1)) +
      geom_point() +
      geom_line()
    
    save_ggplot(plot = pk_plot, title = paste0("Optimal_pK_sample_", sample), folder = parametrization_folder,
                width = 3000, height = 3000)
    
    # Creando el gráfico combinado para de detección de doublets
    
    # DimPlot para visualizar los doublets en UMAP
    plot1 <- DimPlot(seu_sample, group.by = 'doublet_finder', reduction = 'umap', pt.size = 4) + 
      ggtitle("DimPlot Doublet Detection") + theme(plot.title = element_text(hjust = 0.5)) 
    
    # Violin plots de métricas de calidad
    plot2 <- VlnPlot(seu_sample, group.by = 'doublet_finder', features = "nFeature_RNA", pt.size = 0) + 
      ggtitle("nFeature_RNA") + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
    
    plot3 <- VlnPlot(seu_sample, group.by = 'doublet_finder', features = "nCount_RNA", pt.size = 0) + 
      ggtitle("nCount_RNA") + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
    
    plot4 <- VlnPlot(seu_sample, group.by = 'doublet_finder', features = "mitoRatio", pt.size = 0) + 
      ggtitle("percent_mt") + theme(plot.title = element_text(hjust = 0.5))
    
    # Combinar gráficos en una sola figura
    final_plot <- (plot1 / (plot2 | plot3 | plot4)) + 
      patchwork::plot_annotation(title = paste0("Doublet Detection - ", sample))
    
    final_plot <- final_plot + plot_annotation(tag_levels = "A")
    save_ggplot(final_plot, title = paste0("Doublet_Detection_", sample), folder = output_dir,
                width = 4000,
                height = 4000)
    
    doublet_report[i, 3] <- nrow(seu_sample@meta.data[which(seu_sample@meta.data$doublet_finder == 'Singlet'),])
    doublet_report[i, 4] <- nrow(seu_sample@meta.data[which(seu_sample@meta.data$doublet_finder == 'Doublet'),])
    
    print("Removing Doublets from the dataset")
    seu_sample <- subset(seu_sample,
                         subset = doublet_finder == 'Singlet')

    list_seurat[[i]] <- seu_sample
    }# Key for for loop
  
  print("Saving DoubletFinder Report")
  save_dataframe(df = doublet_report, title = "DoubletFinder_Report", folder = output_dir)
  
  
  if(save_intermediates){
    print("Saving Seurat List")
    saveRDS(object = list_seurat, file.path(output_dir, "Seurat_Doublet_Free_Object.rds"))
    
    rm(list_seurat, doublet_wb, doublet_report)
    gc()
    
    list_seurat <- readRDS(file = file.path(output_dir, "Seurat_Doublet_Free_Object.rds"))
  }
  
  print_centered_note("End Of Doublet Detection ")
  
  return(list_seurat)
} # Key to close function
