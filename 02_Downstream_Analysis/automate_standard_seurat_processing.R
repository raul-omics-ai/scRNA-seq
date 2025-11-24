########## 04/02/2025 ##########

#' Automatic Standard Seurat Processing Pipeline
#'
#' @description
#' The `automatic_standard_seurat()` function performs a complete and standardized
#' Seurat-based preprocessing pipeline for single-cell RNA-seq (scRNA-seq) data.
#' It applies normalization, feature selection, scaling, dimensionality reduction,
#' cell cycle scoring (optional), and visualization steps to a list of Seurat objects.
#' Results, plots, and reports are automatically saved in a dedicated directory.
#'
#' @param list_seurat A named list of Seurat objects, each representing one biological sample.
#'   The list should contain objects already pre-filtered and containing basic QC metadata
#'   (e.g. `nFeature_RNA`, `nCount_RNA`, `mitoRatio`, `log10GenesPerUMI`).
#'
#' @param where_to_save Character string. Path to the folder where the output directory
#'   (`Standard_Seurat_Processing_YYYY-MM-DD/`) will be created.  
#'   Default: `getwd()` (the current working directory).
#'
#' @param save_intermediate_files Logical. If `TRUE` (default), saves intermediate objects
#'   (such as the processed Seurat list) to `.rds` files for reproducibility.
#'
#' @param specie Character. Species code to determine whether to perform cell cycle analysis.
#'   Must be `"hsa"` for human data. Default: `"hsa"`.
#'
#' @param cell_cycle_analysis Logical. If `TRUE` (default) and `specie == "hsa"`,
#'   performs cell cycle scoring using predefined S and G2/M gene sets.
#'
#' @details
#' The pipeline consists of several sequential blocks:
#'
#' **1. Loading dependencies and preparing environment**  
#'  - Checks and installs required CRAN packages (`ggplot2`, `Seurat`, `openxlsx`, etc.).  
#'  - Loads custom helper scripts:
#'    - `create_sequential_dir.R`: Creates sequentially numbered output directories.  
#'    - `Automate_Saving_ggplots.R`: Automates saving of ggplot figures.  
#'    - `automate_saving_dataframes_xlsx_format.R`: Saves data frames in Excel format.  
#'    - `print_centered_note_v1.R`: Prints formatted console notes.
#'
#' **2. Standard Seurat preprocessing per sample**  
#'  - Normalization via `NormalizeData()`  
#'  - Highly variable gene selection (`FindVariableFeatures()`)  
#'  - Data scaling (`ScaleData()`)  
#'  - Principal component analysis (`RunPCA()`)
#'  - Automatic PC selection via cumulative variance analysis  
#'  - Elbow plots generated and saved per sample.
#'
#' **3. Source of variation analysis**  
#'  - Adds categorical mitochondrial ratio bins (`percent_mito_CAT()`)  
#'  - Optional: cell cycle phase scoring (S, G2/M) using `CellCycleScoring_edited.R`.
#'
#' **4. Non-linear dimensionality reduction**  
#'  - Automatically selects number of PCs per sample.  
#'  - Runs `RunUMAP()` and `RunTSNE()` accordingly.  
#'  - Generates and saves a “Selected Dimensions Report”.
#'
#' **5. Visualization and reporting**  
#'  - Generates `DimPlot()` visualizations for PCA, UMAP, and tSNE, grouped by `mitoFr`
#'    and optionally `Phase` (cell cycle).  
#'  - Combines plots across samples using `cowplot::plot_grid()` and saves them automatically.  
#'  - Optionally saves the processed list of Seurat objects to disk.
#'
#' @return
#' A list of processed Seurat objects, each containing:
#'  - Normalized data  
#'  - Variable features  
#'  - PCA, UMAP, and tSNE reductions  
#'  - Mitochondrial and (optionally) cell cycle annotations  
#'  
#' Additionally, the function creates a results directory containing:
#'  - `ElbowPlot_*.png` figures  
#'  - UMAP/tSNE visualizations grouped by cell cycle and mitoRatio  
#'  - `Selected_Dimensions_Report.xlsx`  
#'  - `Processed_Seurat_List.rds` (if `save_intermediate_files = TRUE`)
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' processed_list <- automatic_standard_seurat(
#'   list_seurat = list(sample1 = Seurat1, sample2 = Seurat2),
#'   where_to_save = "~/Results/Seurat_Pipeline/",
#'   specie = "hsa",
#'   cell_cycle_analysis = TRUE
#' )
#' }
#'
#' @seealso
#' [Seurat::NormalizeData()], [Seurat::RunPCA()], [Seurat::RunUMAP()],
#' [Seurat::FindVariableFeatures()], [Seurat::ScaleData()], [CellCycleScoring()]
#'
#' @author
#' Raúl — Bioinformatician
#'
#' @export

automatic_standard_seurat <- function(SeuratObject, 
                                      where_to_save = NULL, 
                                      save_intermediate_files = T,
                                      cell_cycle_analysis = T){

  # =================================================================================== #
  # ==== BLOCK 1: LOADING LIBRARIES, CUSTOM FUNCTIONS AND CREATING NEW DIRECTORIES ====
  # =================================================================================== #
  # installing packages 
  print('Instalando los paquetes necesarios para el análisis...')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'dplyr', 'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # creating the saving directory
  where_to_save = ifelse(is.null(where_to_save), getwd(), where_to_save)
  
  saving_directory <- create_sequential_dir(where_to_save, name = "Standard_Seurat_Processing")
    
  # custom functions
  source("~/Documentos/09_scripts_R/create_sequential_dir.R")
  source("~/Documentos/09_scripts_R/Automate_Saving_ggplots.R")
  source("~/Documentos/09_scripts_R/automate_saving_dataframes_xlsx_format.R")
  source("~/Documentos/09_scripts_R/print_centered_note_v1.R")
  
  percent_mito_CAT <- function(seurat_object) {
    # Calcular los cuantiles
    mito_breaks <- quantile(
      seurat_object@meta.data$mitoRatio,
      probs = c(0, 0.25, 0.5, 0.75, 1),
      na.rm = TRUE
    )
    
    # Verificar si los breaks son únicos
    if (length(unique(mito_breaks)) < length(mito_breaks)) {
      # Añadir un pequeño ruido a los cuantiles repetidos para hacerlos únicos
      mito_breaks <- mito_breaks + seq(0, 1e-6, length.out = length(mito_breaks))
    }
    
    # Crear los intervalos usando breaks únicos
    seurat_object$mitoFr <- cut(
      seurat_object@meta.data$mitoRatio,
      breaks = mito_breaks,
      labels = c("Low", "Medium", "Medium high", "High"),
      include.lowest = TRUE
    )
    
    return(seurat_object)
  } # Key of the percentmito function

  # ========================================================== #
  # ==== BLOCK 2: Split Each Biological Sample To Analyze ====
  # ========================================================== #
  print_centered_note("Starting standard seurat processing")
  
  list_seurat <- lapply(list_seurat, function(x){NormalizeData(x)})
  print('NormalizeData OK')
  list_seurat <- lapply(list_seurat, function(x){
    FindVariableFeatures(x,
                         selection.method = "vst",
                         nfeatures = 2000,
                         verbose = FALSE)
    
  })
  print('FindVariableFeatures OK')
  list_seurat <- lapply(list_seurat, function(x){ScaleData(x)})
  print("ScaleData OK")
  list_seurat <- lapply(list_seurat, function(x){RunPCA(x)})
  print("RunPCA OK")
  print("Creating Elbow Plots")
  
  # Calculo automático de las dimensiones a utilizar
  lapply(list_seurat, function(x){
    # Determine percent of variation associated with each PC
    pct <- x[["pca"]]@stdev / sum(x[["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    # Minimum of the two calculation
    pcs <- min(co1, co2)
    
    # Create a dataframe with values
    plot_df <- data.frame(pct = pct, 
                          cumu = cumu, 
                          rank = 1:length(pct))
    
    # Elbow plot to visualize 
    elbowplot <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
      geom_text() + 
      geom_vline(xintercept = 90, color = "grey") + 
      geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
      theme_bw()
    
    save_ggplot(plot = elbowplot, 
                folder = saving_directory,
                title = paste0("ElbowPlot_Sample_", unique(x@meta.data$orig.ident)),
                width = 2000, height = 2000)
    
  })

  # ======================================================= #
  # ==== BLOCK 3: SOURCES OF VARIATION FOR SEURAT LIST ====
  # ======================================================= #
  print_centered_note("Starting analysis unwanted sources of variation")
  # MitoPercent
  print('Percent of mitochondrial damage')
  
  # Aplicar la función a la lista de objetos Seurat
  list_seurat <- lapply(list_seurat, percent_mito_CAT)
  
  # CellCycle
  if(specie == "hsa" & cell_cycle_analysis){
      print('Starting Cell Cycle analysis')
      # cell cycle markers
      
      s_genes <- c("UBR7",     "RFC2",     "RAD51",    "MCM2",     "TIPIN",    "MCM6",     
                   "UNG",      "POLD3",    "WDR76",    "CLSPN",    "CDC45",   "CDC6",
                   "MSH2",     "MCM5",     "POLA1",    "MCM4",     "RAD51AP1", "GMNN",
                   "RPA2",     "CASP8AP2", "HELLS",    "E2F8",    "GINS2",    "PCNA",
                   "NASP",     "BRIP1",    "DSCC1",    "DTL",     "CDCA7",    "CENPU",
                   "ATAD2",    "CHAF1B",   "USP1",    "SLBP",     "RRM1",     "FEN1",
                   "RRM2",     "EXO1",     "CCNE2",   "TYMS",     "BLM",      "PRIM1",
                   "UHRF1")
      
      g2m_genes <- c("NCAPD2",  "ANLN",    "TACC3",   "HMMR",    "GTSE1",   "NDC80",
                     "AURKA",   "TPX2",    "BIRC5",   "G2E3",    "CBX5",   "RANGAP1",
                     "CTCF",    "CDCA3",   "TTK",     "SMC4",    "ECT2",    "CENPA",
                     "CDC20",   "NEK2",    "CENPF",   "TMPO",    "HJURP",   "CKS2",   
                     "DLGAP5",  "PIMREG",  "TOP2A",   "PSRC1",   "CDCA8",   "CKAP2",
                     "NUSAP1",  "KIF23",   "KIF11",   "KIF20B",  "CENPE",   "GAS2L3", 
                     "KIF2C",   "NUF2",    "ANP32E",  "LBR",     "MKI67",   "CCNB2",
                     "CDC25C",  "HMGB2",   "CKAP2L",  "BUB1",    "CDK1",    "CKS1B",  
                     "UBE2C",   "CKAP5",   "AURKB",   "CDCA2",   "TUBB4B",  "JPT1")
      
      source("~/Documentos/09_scripts_R/CellCycleScoring_edited.R")
      list_seurat <- lapply(list_seurat, function(x){
        CellCycleScoring_edited(x, 
                                g2m.features = g2m_genes, 
                                s.features = s_genes)
      })
      print('CellCycleScoring OK')
    } # if cell cycle analysis key

  # ====================================================================== #
  # ==== BLOCK 4: NON-LINEAR DIMENSIONALITY REDUCTION FOR SEURAT LIST ====
  # ====================================================================== #
  print_centered_note('Starting with Non-Linear Dimensionality Reduction Step')
  
  save_dims_per_sample <- c()
  dimensions_report <- data.frame(matrix(nrow=0,ncol=2))
  colnames(dimensions_report) = c('Sample', 'Number of PC')
  
  for(sample in 1:length(list_seurat)){
    print(paste0("Computing PCs for sample ", names(list_seurat)[sample]))
    # Calculo automático de las dimensiones a utilizar
    pct <- list_seurat[[sample]][["pca"]]@stdev / sum(list_seurat[[sample]][["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    # Minimum of the two calculation
    pcs <- min(co1, co2)
    dims_per_sample <- pcs
    save_dims_per_sample <- c(save_dims_per_sample, dims_per_sample)
    print(paste0('A total of ', pcs, ' PCs were selected for UMAP and TSNE'))
    list_seurat[[sample]] <- RunUMAP(list_seurat[[sample]], dims = 1:dims_per_sample)
    print('RunUMAP OK')
    list_seurat[[sample]] <- RunTSNE(list_seurat[[sample]], dims = 1:dims_per_sample)
    print('RunTSNE OK')
    
    dimensions_report[sample,1] = names(list_seurat)[sample]
    dimensions_report[sample,2] = dims_per_sample
  } # for loop for computing PCs, UMAP and plots key
  
  save_dataframe(df = dimensions_report, 
                 title = "Selected_Dimensions_Report", 
                 folder = saving_directory)

  # ================================ #
  # === BLOCK 5: VISUALIZATIONS ====
  # ================================ #
  print_centered_note("Creating Visualizations About Sources of Variation")
  num_samples <- length(list_seurat)
  
  i <- 1
  groups <- ifelse(cell_cycle_analysis & specie == "hsa", c("Phase", "mitoFr"), "mitoFr")
  for (reduction in c("pca", "umap", "tsne")) {
    for (group in groups) {
      plotlist <- lapply(list_seurat, function(x) {
        DimPlot(x, reduction = reduction, group.by = group) + ggtitle(x$orig.ident)
      })
      
      # Determinar número de columnas dinámicamente
      ncol_value <- if (num_samples <= 3) {
        num_samples
      } else if (num_samples %% 4 == 0) {
        2
      } else {
        3
      }
      
      plot <- plot_grid(plotlist = plotlist, ncol = ncol_value)
      save_ggplot(plot, title = paste0(reduction, "_", group, "_All_Samples"), 
                  folder = saving_directory,
                  dpi = 300, width = 3000, height = 2000)
      
      i <- i + 1 # Incrementar el contador
    } # for loop over groups key
  } # for loop over reductions key
  
  
  if(save_intermediate_files){
    print("Saving Seurat List")
    saveRDS(list_seurat, file = file.path(saving_directory, "Processed_Seurat_List.rds"))
    
    rm(list_seurat, i, num_samples, dimensions_report)
    gc()
    list_seurat <- readRDS(file = file.path(saving_directory, "Processed_Seurat_List.rds"))
  }
  
  print_centered_note("End Of The Standard Seurat Processing Pipeline ")
  
  return(list_seurat)
  
  } # Key of the function
