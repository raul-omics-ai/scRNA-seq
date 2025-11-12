########## 7/5/2024 ##########
#------------------------------------------------------------------------------#
########### ANÁLISIS DE CONTROL DE CALIDAD PARA DATOS DE SCRNA-SEQ #############
#------------------------------------------------------------------------------#

# Este script lo he creado para obtener las métricas del control de calidad de los
# datos de scRNA-seq. 

# INPUT: Una lista de objetos Seurat (cada objeto dentro de la lista será una
# muestra)

# OUTPUT: Se devolverá la lista de objetos seurat con los metadatos nuevos y
# una carpeta que constendrá los siguientes archivos:

#   - Un archivo metadata.xlsx con la información de los metadatos recogida de todas
#     las células
#   - Un archivo plot_interpretation.xlsx con información para interpetar los
#     gráficos
#   - Distintos gráficos sobre la calidad de los datos de cada muestra de scRNA-seq


automatic_qc_scrnaseq <- function(list_srn, 
                                  where_to_save = NULL,
                                  specie = "hsa"){
  print('Comenzando el análisis de control de calidad de los datos de scRNA-seq...')
  
  ########## FUNCIONES ##########
  if(specie == "hsa"){
    pat <- "^MT-"
  }
  if(specie == "mmu"){
    pat = "^mt-"
  }
  
  percent_mito <- function(seurat_object){
    seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, pattern = pat)
    seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100
    return(seurat_object)
  }
  
  ########## INSTALLING PACKAGES ##########
  list.of.packages = c('ggplot2', 'Seurat', 'SeuratObject','openxlsx')
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  start.time <- Sys.time()
  
  
  # Diretorio donde se van a guardar todas las imágenes
  if(is.null(where_to_save)){
    where_to_save <- './'
  }
  
  results_folder <- paste0(where_to_save, '1_QUALITY_CONTROL_', Sys.Date())
  dir.create(results_folder)
  setwd(results_folder)
  # Vamos a ver si todas las columnas de los metadatos tienen el mismo número de columnas
  num_columns_list <- (lapply(list_srn, function(x){num_columns <- ncol(x@meta.data)}))
  if (length(unique(num_columns_list)) == 1) {
    print("Todos los objetos tienen el mismo número de columnas en los metadatos.")
    print(paste("Número de columnas:", num_columns_list[[1]]))
  }
  
  # Compute mitoRatio column (values from 0 to 1) to each seurat_object
  list_srn <- lapply(list_srn, percent_mito)
  
  # Create the combined_metadata dataframe to perform qc metrics
  print('Creando los metadatos...')
  combined_metadata <- data.frame(matrix(nrow = 0, ncol = num_columns_list[[1]]))
  for(sample in names(list_srn)){
    combined_metadata <- rbind(combined_metadata, list_srn[[sample]]@meta.data)
    
  }
  
  combined_metadata$log10GenesPerUMI <- log10(combined_metadata$nFeature_RNA) / log10(combined_metadata$nCount_RNA)
  
  for(i in names(list_srn)){
    list_srn[[i]]@meta.data <- combined_metadata[which(combined_metadata$orig.ident == i),]
  }
  print('Guardando los metadatos...')
  write.xlsx(combined_metadata, file=paste0("./metadata.xlsx"), overwrite = T)
  
  ########### PLOTS ##########
  print('Creando el archivo de interpretar los gráficos...')
  interpretation <- data.frame(matrix(nrow=7, ncol=2))
  colnames(interpretation) <- c('QC_Metric', 'Interpretation')
  interpretation$QC_Metric <- c("Células por muestra", 'nUMIs por célula', 'nGenes detectados por célula', 'mitoRatio', 'Complejidad', 'UMIs vs Genes detectados', 'Referencia')
  
  ## Células por muestra
  interpretation[1,2] <- 'Representa el número de células que hay por muestra. Si las muestras son homogéneas, estarán distribuidas iguales en el density plot'
  ## nUMIs por célula
  interpretation[2,2] <- 'nUMIs (nFeatures) hace referencia al número de tránscritos (features) que se han detectado en la secuenciación. Al menos tiene que haber más de 500 UMIs por célula. Si es superior, esto nos indica que es bueno y que podemos utilizarlo porque la secuenciación ha tenido mayor profundidad'
  ## Genes detectados por célula
  interpretation[3,2] <- 'nGene (es lo mismo que nCount) en el número de genes detectados por célula. Se representa de forma también de density plot en el cual lo ideal es que hubiera únicamente un único pico. Si hay un hombro en el lado derecho o una distribución bimodal puede deberse a dos cosas: 1) Que las células han muerto por alguna razón 2) There are biologically different types of cells (poblaciones quiescentes, poco complejas), es decir, puede haber un tipo celular más pequeño que otro (por ejemplo células de gran tamaño pueden tener mayores counts)'
  ## MitoRatio
  interpretation[4,2] <- 'mitoRatio es la medida de la muerte celular expresada de 0 a 1 por la secuenciación de los genes mitocondriales. Dentor de este ratio, los valores aceptables van del 5 al 20 % (0.05 - 0.2)'
  ##Complejidad
  interpretation[5,2] <- 'La complejidad en este contexto hace referencia a la profunidad de secuenciación, en donde las células peor secuenciadas tienen una complejidad mayor debido a que no se ha saturado la secuenciación en ningún gen de las muestras. Las células outliers podrían ser células con una menor complejidad de RNA que otras células y a veces, se pueden detecta contaminacones de otros tipos celulares con baja complejidad (como por ejemplo glóbulos rojos). Generalmente, se espera un valor superior a 0.3 y para eso hemos calculado el log10 de genes por UMI.
  Se calcula como log10 del número de genes entre el log10 del número de tránscritos detectados (UMIs)'
  ##UMIs vs Genes detectados
  interpretation[6,2] <- 'Correlación entre el número de UMIs (tránscritos) y genes detectados. Estas gráficas nos van a permitir crear un filtro para eliminar las células de baja calidad, porque serán las que tengan muy pocos genes y pocos UMIs por célula (se corresponden con el cuadrado inferior izquierdo de la gráfica). Células de buena calidad tendrán un gran número de genes y un gran número de tránscritos (UMIs).

Esta gráfica además también nos permite distinguir aquellas células que tienen un alto número de UMIs (tránscritos) pero un bajo número de genes. Estas células pueden ser también ceúlas que se están muriendo, pero pueden represnetar una oblación de tipos celulares de baja complejidad (como los glóbulos rojos)'
  ## Referencia
  interpretation[7,2] <- 'https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html'
  print('Guardando el archivo de interpretación de los gráficos...')
  write.xlsx(interpretation, file = paste0(where_to_save, results_folder, 'plots_interpretation.xlsx'), overwrite = T, asTable = T)
  print('Creando los gráficos de calidad...')
  ########## CELLS PER SAMPLE ##########
  print('1...Células por muestra')
  png(filename = paste0("./number_cells_per_sample.png"), res = 300, width = 2500, height = 3000)
  print(combined_metadata %>%
          ggplot(aes(x = orig.ident, fill = orig.ident))+
          geom_bar()+
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells"))
  dev.off()
  
  ########## UMI COUNTS PER CELL ##########
  print('2...UMIs por célula')
  png(filename = paste0("./density_plot_UMI_counts_per_cell_per_sample.png"), res = 300, width = 2500,  height = 3000)
  print(combined_metadata %>% 
          ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          scale_x_log10() + 
          theme_classic() +
          ylab("Cell density") +
          xlab("nUMI") +
          geom_vline(xintercept = 500)+
          facet_grid(orig.ident ~.))
  dev.off()
  
  ########## GENES DETECTED PER CELL ##########
  print('3...Genes detectados por célula (density plot)')
  png(paste0("./density_plot_genes_detected_per_cell_per_sample.png"), res = 300, width = 3000, height = 3000)
  print(combined_metadata %>% 
          ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          xlab("nGene") +
          geom_vline(xintercept = 300)+
          facet_grid(orig.ident ~ .))
  dev.off()
  
  # Visualize the distribution of genes detected per cell via boxplot
  print('4...Genes detectados por célula (boxplot)')
  png(paste0("./boxplot_genes_detected_per_cell.png"), res = 300, width = 3000, height = 3000)
  print(combined_metadata %>% 
          ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
          geom_boxplot() + 
          theme_classic() +
          ylab("nGene")+
          #geom_signif(comparisons = list(c("APG1", "APG2"),
          #                               c('APG2', 'APG3'),
          #                               c('APG1', 'APG3')),
          #y_position = c(3.5, 3.5, 3.75))+
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells vs NGenes"))
  dev.off()
  
  ########## UMIS vs GENES DETECTED ##########
  print('5...UMIs vs Genes detectados')
  png(paste0("./scatterplot_genes_vs_umis_per_sample.png"), res = 300, width = 3000, height = 3000)
  print(combined_metadata %>% 
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
          facet_wrap(~orig.ident))
  dev.off()
  
  ########## MITOCONDRIAL COUNTS RATIO ##########
  print('6...Mitocondrial ratio')
  png(paste0("./density_plot_mitochondrial_counts.png"), res = 300, width = 3000, height = 3000)
  print(combined_metadata %>% 
          ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
          geom_density(alpha = 0.2) + 
          scale_x_log10() + 
          theme_classic() +
          geom_vline(xintercept = 0.2)+
          geom_vline(xintercept = 0.05)+
          facet_grid(orig.ident ~ .))
  dev.off()
  
  ########## COMPLEXITY ##########
  print('7...Complejidad')
  png(paste0("./density_plot_Complexity.png"), res = 300, width = 3000, height = 3000)
  print(combined_metadata %>%
          ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
          geom_density(alpha = 0.2) +
          theme_classic() +
          geom_vline(xintercept = 0.8)+
          facet_grid(orig.ident ~.))
  dev.off()
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  setwd("../")
  print('Script terminado!')
  print(paste0('El script ha tardado ', time.taken, ' segundos'))
  print('Nota: Recuerda que los  metadatos se han modificado en cada uno de los objetos seurat.')
  return(list_srn)
  
}



