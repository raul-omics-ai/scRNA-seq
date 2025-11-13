########## 7/5/2024 ##########
#' Automatic Quality Control for scRNA-seq Data
#'
#' @description
#' This function performs an automated **Quality Control (QC)** analysis for
#' **single-cell RNA sequencing (scRNA-seq)** data from a list of `Seurat` objects.
#'
#' It computes key QC metrics such as the number of detected genes, UMIs,
#' mitochondrial gene ratio, and transcriptome complexity.
#' The function also generates several diagnostic plots and exports
#' all results (metadata tables and figures) to structured output folders.
#'
#' @param list_srn 
#' A list of `Seurat` objects.  
#' Each element of the list represents an independent scRNA-seq sample.
#'
#' @param where_to_save 
#' Path to the directory where all output files and figures will be saved.  
#' If `NULL`, results will be saved in the current working directory (`"."`).
#'
#' @param specie 
#' The organism from which the data originates.  
#' Used to identify mitochondrial genes.  
#' Accepted values:
#' \itemize{
#'   \item `"hsa"` for *Homo sapiens* (pattern `"^MT-"`)
#'   \item `"mmu"` for *Mus musculus* (pattern `"^mt-"`)
#' }
#' Default: `"hsa"`.
#'
#' @details
#' The function executes the following main steps:
#'
#' **1. Environment Setup**
#' - Checks for required R packages (`ggplot2`, `Seurat`, `openxlsx`, etc.) and installs any missing ones.
#' - Loads external custom utility scripts for text printing, directory creation, and saving outputs.
#'
#' **2. Preprocessing and Metadata Generation**
#' - Computes the mitochondrial gene ratio (`mitoRatio`) for each cell.
#' - Combines metadata from all samples into a single data frame (`combined_metadata`).
#' - Calculates the `log10GenesPerUMI` metric as a measure of sequencing complexity.
#'
#' **3. Exporting Metadata and Interpretation Notes**
#' - Exports the combined metadata to `Combined_Metadata.xlsx`.
#' - Creates an interpretation table (`Plot_Interpretation_Notes.xlsx`) describing each QC metric and its biological meaning.
#'
#' **4. Generating QC Visualizations**
#' The function creates and saves individual plots for the following metrics:
#' - Number of cells per sample.
#' - Distribution of UMIs per cell.
#' - Distribution and boxplot of detected genes per cell.
#' - Correlation between UMIs and detected genes.
#' - Mitochondrial ratio distribution.
#' - Transcriptome complexity (log10GenesPerUMI).
#'
#' Additionally, a combined figure (`QC_Figures.png`) containing all main plots is created.
#'
#' **5. Output**
#' - Returns the original list of `Seurat` objects with updated metadata containing:
#'   - `mitoRatio`: mitochondrial gene proportion (0–1)
#'   - `log10GenesPerUMI`: complexity metric
#'
#' - Creates an output directory named:
#'   `"Quality_Control"`, containing:
#'   \itemize{
#'     \item `Combined_Metadata.xlsx` — combined metadata table
#'     \item `Plot_Interpretation_Notes.xlsx` — interpretation guide for the QC plots
#'     \item Folder `Individual_Figures/` — individual QC plots
#'     \item `QC_Figures.png` — combined figure with all QC panels
#'   }
#'
#' @return
#' A list of `Seurat` objects updated with additional metadata columns:
#' \itemize{
#'   \item `mitoRatio`: proportion of mitochondrial genes
#'   \item `log10GenesPerUMI`: sequencing complexity metric
#' }
#'
#' The function also saves metadata tables and QC plots in the specified directory.
#'
#' @dependencies
#' Required packages:
#' - `ggplot2`, `Seurat`, `SeuratObject`, `openxlsx`, `patchwork`, `dplyr`
#'
#' External scripts required:
#' - `print_centered_note_v1.R`
#' - `create_sequential_dir.R`
#' - `automate_saving_dataframes_xlsx_format.R`
#' - `Automate_Saving_ggplots.R`
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' qc_results <- automatic_qc_scrnaseq(
#'   list_srn = list(sample1 = seurat_obj1, sample2 = seurat_obj2),
#'   where_to_save = "~/Results/QC_scRNAseq/",
#'   specie = "hsa"
#' )
#'
#' # Access updated metadata
#' head(qc_results[[1]]@meta.data)
#' }
#'
#' @references
#' HBC Training:  
#' <https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html>
#'
#' @author
#' Raúl (Bioinformatician)
#'
#' @export

automatic_qc_scrnaseq <- function(list_srn, 
                                  where_to_save = NULL,
                                  specie = "hsa"){
  ########## INSTALLING PACKAGES ##########
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx')
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  start.time <- Sys.time()

  ########## CUSTOM FUNCTIONS ##########
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")

  percent_mito <- function(seurat_object, specie){
    if(specie == "hsa"){
      pat <- "^MT-"
    }
    if(specie == "mmu"){
      pat <- "^mt-"
    }
    seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, pattern = pat)
    seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100
    return(seurat_object)
  }

  # ============================================== #
  # ==== STEP 1: SETTING UP WORKING DIRECTORY ====
  # ============================================== #
  print_centered_note('Starting the QC of scRNA-seq data')
  
  # Working directory
  where_to_save <- ifelse(is.null(where_to_save), ".", where_to_save)
  results_folder <- create_sequential_dir(path = where_to_save, name = "Quality_Control")

  
  # Checkpoint: Check if all samples have the same number of columns in metadata
  num_columns_list <- (lapply(list_srn, function(x){num_columns <- ncol(x@meta.data)}))
  if (length(unique(num_columns_list)) == 1) {
    print("All objects have the same number of columns in metadata")
    print(paste("Number of columns:", num_columns_list[[1]]))
  }
  
  # Compute mitoRatio column (values from 0 to 1) to each seurat_object
  list_srn <- lapply(list_srn, percent_mito)
  
  # Create the combined_metadata dataframe to perform qc metrics
  print('Combining metadata')
  combined_metadata <- data.frame(matrix(nrow = 0, ncol = num_columns_list[[1]]))
  for(sample in names(list_srn)){
    combined_metadata <- rbind(combined_metadata, list_srn[[sample]]@meta.data)
    
  }
  
  combined_metadata$log10GenesPerUMI <- log10(combined_metadata$nFeature_RNA) / log10(combined_metadata$nCount_RNA)
  
  for(i in names(list_srn)){
    list_srn[[i]]@meta.data <- combined_metadata[which(combined_metadata$orig.ident == i),]
  }
  print('Saving metadata')
  save_dataframe(combined_metadata, title = "Combined_Metadata", folder = results_folder)

  # ===================================== #
  # ==== STEP 2: INTERPRETATION FILE ====
  # ===================================== #
  print_centered_note('Creating the file to interpret plots')
  interpretation <- data.frame(matrix(nrow=7, ncol=2))
  colnames(interpretation) <- c('QC_Metric', 'Interpretation')
  interpretation$QC_Metric <- c("Cells per sample", 'UMIs per cell', 'nGenes detected per cell', 
                                'mitoRatio', 'Complexity', 'UMIs vs Genes detected', 'Reference')
  
  ## Cells per sample
  interpretation[1,2] <- 'Represents the number of cells per sample. If the samples are homogeneous, they will be distributed evenly in the density plot.'
  ## nUMIs per cell
  interpretation[2,2] <- 'nUMIs (same as nFeatures) refers to the number of transcripts (features) that have been detected in the sequencing. There must be at least 500 UMIs per cell. If there are more, this indicates that it is good and that we can use it because the sequencing has had greater depth.'
  ## Genes detected per cell
  interpretation[3,2] <- 'nGene (same as nCount) is the number of genes detected per cell. It is also represented as a density plot, which ideally should have only a single peak. If there is a shoulder on the right side or a bimodal distribution, this may be due to two things: 1) The cells have died for some reason. 2) There are biologically different types of cells (quiescent, less complex populations), i.e., there may be one cell type that is smaller than another (for example, large cells may have higher counts).'
  ## MitoRatio
  interpretation[4,2] <- 'mitoRatio is a measure of cell death expressed from 0 to 1 by sequencing mitochondrial genes. Within this ratio, acceptable values range from 5 to 20% (0.05 - 0.2).'
  ## Complexity
  interpretation[5,2] <- 'Complexity in this context refers to sequencing depth, where poorly sequenced cells have greater complexity because sequencing has not been saturated in any gene in the samples. Outlier cells could be cells with lower RNA complexity than other cells, and sometimes contamination from other cell types with low complexity (such as red blood cells) can be detected. Generally, a value greater than 0.3 is expected, and for this we have calculated the log10 of genes per UMI.
  It is calculated as log10 of the number of genes divided by log10 of the number of detected transcripts (UMIs).'
  ##UMIs vs Genes
  interpretation[6,2] <- 'Correlation between the number of UMIs (transcripts) and detected genes. These graphs will allow us to create a filter to eliminate low-quality cells, because they will be the ones with very few genes and few UMIs per cell (corresponding to the lower left square of the graph). Good-quality cells will have a large number of genes and a large number of transcripts (UMIs).

This graph also allows us to distinguish those cells that have a high number of UMIs (transcripts) but a low number of genes. These cells may also be dying cells, but they may represent an obligation of low-complexity cell types (such as red blood cells).'
  ## Reference
  interpretation[7,2] <- 'https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html'
  
  print('Saving the graphics interpretation file')
  
  save_dataframe(df = interpretation, title = "Plot_Interpretation_Notes", folder = results_folder)

  # ========================================= #
  # ==== STEP 3: CREATING VISUALIZATIONS ====
  # ========================================= #
  print_centered_note('Creating QC plots')
  individual_images_dir <- create_sequential_dir(path = results_folder, name = "Individual_Figures")

  ########## CELLS PER SAMPLE ##########
  cat('\n1...Cells per sample\n')
  
  cells_per_sample <- combined_metadata %>%
    ggplot(aes(x = orig.ident, fill = orig.ident))+
    geom_bar()+
    xlab(NULL)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("Number Cells Per Sample")+
    theme(legend.position = "none")
  
  save_ggplot(plot = cells_per_sample, title = "Number_Cells_Per_Sample", 
              folder = individual_images_dir,
              width = 2500, height = 3000)
  
  ########## UMI COUNTS PER CELL ##########
  cat('\n2...UMIs per cell\n')
  umi_per_cell <- combined_metadata %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    xlab("nUMI") +
    geom_vline(xintercept = 500) +
    facet_grid(orig.ident ~.) + 
    ggtitle("UMIs per Cell by Sample")+
    theme(legend.position = "none")
  
  
  save_ggplot(plot = umi_per_cell, title = "Density_Plot_UMI_per_Cell_Per_Sample",
              folder = individual_images_dir, width = 2500, height = 3000)
    
  ########## GENES DETECTED PER CELL ##########
  cat('\n3...Genes detected per cell (density plot)\n')
  Genes_per_cell <- combined_metadata %>% 
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    xlab("nGene") +
    geom_vline(xintercept = 300)+
    facet_grid(orig.ident ~ .)+
    ggtitle("Genes Detected per Cell by Sample")+
    theme(legend.position = "none")
  
  save_ggplot(plot = Genes_per_cell, title = "Density_Plot_Genes_Detected_per_Cell_Per_Sample",
              folder = individual_images_dir,
              width = 3000, height = 3000)
  
  # Visualize the distribution of genes detected per cell via boxplot
  cat('\n4...Genes detected per cell (boxplot)\n')
  boxplot_genes_per_cell <- combined_metadata %>% 
    ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
    geom_boxplot() + 
    theme_classic() +
    ylab("nGene")+
    xlab(NULL)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells vs NGenes")+
    theme(legend.position = "none")
  
  save_ggplot(plot = boxplot_genes_per_cell, title = "Boxplot_Genes_Detected_Per_Cell",
                 folder = individual_images_dir, width = 3000, height = 3000)
    
  ########## UMIS vs GENES DETECTED ##########
  cat('\n5...UMIs vs Genes detected\n')
  umis_vs_genes <- combined_metadata %>% 
                     ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
                     geom_point() + 
                     scale_colour_gradient(low = "gray90", high = "black") +
                     stat_smooth(method=lm) +
                     scale_x_log10() + 
                     scale_y_log10() + 
                     theme_classic() +
                     xlab("nUMI")+
                     ylab("nGene")+
                     geom_vline(xintercept = 500) +
                     geom_hline(yintercept = 250) +
                     facet_wrap(~orig.ident)+
                     ggtitle("UMI vs Genes Correlation")
  
  
  save_ggplot(plot = umis_vs_genes, title = "Scatterplot_Genes_VS_UMIs", folder = individual_images_dir,
              width = 3000, height = 3000)
  
  ########## MITOCONDRIAL COUNTS RATIO ##########
  cat('\n6...Mitocondrial ratio\n')
  mito_ratio_plot <- combined_metadata %>% 
    ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)+
    geom_vline(xintercept = 0.05)+
    facet_grid(orig.ident ~ .)+
    ggtitle("MitoRatio")+
    theme(legend.position = "none")
  
  save_ggplot(plot = mito_ratio_plot, title = "Density_Plot_MitoRatio", folder = individual_images_dir,
              width = 3000, height = 3000)
  
  ########## COMPLEXITY ##########
  cat('\n7...Complexity\n')
  complexity <- combined_metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)+
    facet_grid(orig.ident ~.)+
    ggtitle("Complexity")+
    theme(legend.position = "bottom")
  
  save_ggplot(plot = complexity, title = "Density_Plot_Complexity", folder = individual_images_dir,
              width = 3000, height = 3000)

  ########## COMBINED PLOT ##########
  cat('\n8...Combined Figure\n')
  stack1 <- (cells_per_sample / boxplot_genes_per_cell) | 
    (umi_per_cell / Genes_per_cell / mito_ratio_plot / complexity) | umis_vs_genes
  
  stack1 <- stack1 + plot_annotation(tag_levels = 'A', 
                                    theme = theme(plot.tag = element_text(face = "bold")))
  
  save_ggplot(plot = stack1, title = "QC_Figures", folder = results_folder, width = 5000,
              height = 4000)
 
  
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)

  print_centered_note('End of the script')
  print(paste0('This script takes ', time.taken, ' seconds'))
  print('Note:  Remember that the metadata has been modified in each of the seurat objects.')
  return(list_srn)

}






