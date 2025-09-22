########## 04/02/2025 ##########

# ============================================ #
# AUTOMATIC STANDAR SEURAT PROCESSING PIPELINE #
# ============================================ #

# This sript will be used for create a new function to perform the standar seurat pipeline for processing
# scRNA-seq data.

automatic_standard_seurat <- function(SeuratObject, 
                                      where_to_save = NULL, 
                                      save_intermediate_files = T,
                                      cell_cycle_analysis = T){
  # BLOCK 1: LOADING LIBRARIES, CUSTOM FUNCTIONS AND CREATING NEW DIRECTORIES ####
  # installing packages 
  print('Instalando los paquetes necesarios para el análisis...')
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx', 'parallel', 
                       'DoubletFinder', 'harmony', 'dplyr', 'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # creating the saving directory
  if(is.null(where_to_save)){
    where_to_save <- getwd()
  }
  
  saving_directory <- file.path(where_to_save, "03_STANDARD_PROCESSING")
  dir.create(saving_directory, recursive = T, showWarnings = F)
  
  # custom functions
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
  
  # header
  separator <- paste0(rep("-", 80), collapse = "")
  cat("\n")
  print(separator)
  print("                 STARTING WITH THE STANDARD PIPELINE OF SEURAT                  ")
  print(separator)
  
  # BLOCK 2: SPLIT EACH  BIOLOGICAL SAMPLE TO ANALYZE ####
  print("Splitting Seurat Object to analyze each biological sample individually")
  SeuratList <- SplitObject(SeuratObject, split.by = "orig.ident")
  
  print("Starting the standard pipeline of Seurat for each biological sample")
  SeuratList <- lapply(SeuratList, function(x){NormalizeData(x)})
  print('NormalizeData OK')
  SeuratList <- lapply(SeuratList, function(x){
    FindVariableFeatures(x,
                         selection.method = "vst",
                         nfeatures = 2000,
                         verbose = FALSE)
    
  })
  print('FindVariableFeatures OK')
  SeuratList <- lapply(SeuratList, function(x){ScaleData(x)})
  print("ScaleData OK")
  SeuratList <- lapply(SeuratList, function(x){RunPCA(x)})
  print("RunPCA OK")
  print("Calculando los elbow plots...")
  
  # Calculo automático de las dimensiones a utilizar
  lapply(SeuratList, function(x){
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
    png(filename = file.path(saving_directory, paste0("Elbowplot_Sample_", unique(x@meta.data$orig.ident), ".png")), res = 300, width = 2000, height = 2000)
    print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
            geom_text() + 
            geom_vline(xintercept = 90, color = "grey") + 
            geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
            theme_bw())
    dev.off()
  })
  
  # BLOCK 3: SOURCES OF VARIATION FOR SEURAT LIST ####
  # MitoPercent
  print('Starting analysis for percent mito')
  
  # Aplicar la función a la lista de objetos Seurat
  SeuratList <- lapply(SeuratList, percent_mito_CAT)
  
  # CellCycle
  if(cell_cycle_analysis){
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
      SeuratList <- lapply(SeuratList, function(x){
        CellCycleScoring_edited(x, 
                                g2m.features = g2m_genes, 
                                s.features = s_genes)
      })
      print('CellCycleScoring OK')
    } # if cell cycle analysis key
  
  # BLOCK 4: NON-LINEAR DIMENSIONALITY REDUCTION FOR SEURAT LIST ####
  print('Starting non-linear dimensionality reduction for seurat')
  save_dims_per_sample <- c()
  dimensions_report <- data.frame(matrix(nrow=0,ncol=2))
  colnames(dimensions_report) = c('Sample', 'Number of PC')
  
  for(sample in 1:length(SeuratList)){
    print(paste0("Computing PCs for sample ", names(SeuratList)[sample]))
    # Calculo automático de las dimensiones a utilizar
    pct <- SeuratList[[sample]][["pca"]]@stdev / sum(SeuratList[[sample]][["pca"]]@stdev) * 100
    # Calculate cumulative percents for each PC
    cumu <- cumsum(pct)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co1 <- which(cumu > 90 & pct < 5)[1]
    co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
    # Minimum of the two calculation
    pcs <- min(co1, co2)
    dims_per_sample <- pcs
    save_dims_per_sample <- c(save_dims_per_sample, dims_per_sample)
    print(paste0('A total of  ', pcs, ' were selected for UMAP and TSNE'))
    SeuratList[[sample]] <- RunUMAP(SeuratList[[sample]], dims = 1:dims_per_sample)
    print('RunUMAP OK')
    SeuratList[[sample]] <- RunTSNE(SeuratList[[sample]], dims = 1:dims_per_sample)
    print('RunTSNE OK')
    
    dimensions_report[sample,1] = names(SeuratList)[sample]
    dimensions_report[sample,2] = dims_per_sample
  } # for loop for computing PCs, UMAP and plots key
  write.xlsx(dimensions_report, file = file.path(saving_directory, "Selected_Dimensions_Report.xlsx"), overwrite = T)
  
  # BLOCK 5: VISUALIZATIONS
  print("Creating Visualizations about Sources of Variation")
  num_samples <- length(SeuratList)
  
  i <- 1
  for (reduction in c("pca", "umap", "tsne")) {
    for (group in c("Phase", "mitoFr")) {
      title <- sprintf("%02d_%s_%s_All_Samples.png", i, toupper(reduction), group)
      plotlist <- lapply(SeuratList, function(x) {
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
      
      png(filename = file.path(saving_directory, title), res = 300, width = 2000, height = 2000)
      print(plot_grid(plotlist = plotlist, ncol = ncol_value))
      dev.off()
      
      i <- i + 1 # Incrementar el contador
    } # for loop over groups key
  } # for loop over reductions key
  
  print("Merging Seurat List")
  SeuratObject <- merge(x = SeuratList[[1]],
                        y = SeuratList[2:length(SeuratList)],
                        merge.data = TRUE, add.cell.ids = names(SeuratList))
  
  
  if(save_intermediate_files){
    print("Saving Seurat List")
    saveRDS(SeuratObject, file = file.path(saving_directory, "Processed_Seurat_List.rds"))
    
    rm(SeuratObject, i, num_samples, dimensions_report)
    gc()
    SeuratObject <- readRDS(file = file.path(saving_directory, "Processed_Seurat_List.rds"))
  }
  
  
  print(separator)
  print("                 END OF THE STANDAR SEURAT PROCESSING PIPELINE                  ")
  print(separator)
  
  return(SeuratObject)
  
  } # Key of the function
