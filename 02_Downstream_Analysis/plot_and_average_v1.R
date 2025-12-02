########## 02/12/2025 ##########

#' Plot Marker Expression and Compute Average Expression in a Seurat Object
#'
#' This function generates feature plots for a set of marker genes using a Seurat
#' object and computes the average expression of those markers across all cell groups.
#' Results are returned as a data frame containing the average expression per gene
#' and the overall mean expression across groups.
#'
#' @param seurat_obj A \code{Seurat} object containing the data to be visualized
#'   and analyzed.
#' @param markers A character vector specifying the marker genes to plot and
#'   include in the average expression calculation.
#' @param palette A color palette passed to \code{FeaturePlot}. Defaults to
#'   \code{magma(10)}.
#' @param min_cutoff Minimum expression cutoff for \code{FeaturePlot}. Defaults
#'   to \code{"q10"}.
#' @param max_cutoff Maximum expression cutoff for \code{FeaturePlot}. Defaults
#'   to \code{"q90"}.
#' @param assay Name of the assay from which to compute average expression.
#'   Defaults to \code{"RNA"}.
#' @param reduction Dimensionality reduction to use for plotting. Defaults to
#'   \code{"umap"}.
#'
#' @details
#' The function first ensures that all required packages are installed and loaded.
#' It then generates feature plots for the specified markers using Seurat's
#' \code{FeaturePlot}. The number of columns in the plot grid is automatically
#' determined based on the number of markers.
#'
#' Next, the function computes average expression using \code{AverageExpression}
#' and formats the output into a data frame. A new column, \code{media}, is added
#' containing the mean expression across all groups for each marker.
#'
#' @return A data frame with one row per gene, including the average expression
#'   across cell groups and an additional column \code{media} indicating the mean
#'   expression across all groups.
#'
#' @examples
#' \dontrun{
#' markers <- c("GeneA", "GeneB")
#' avg_df <- plot_and_average(seurat_obj = my_seurat,
#'                            markers = markers,
#'                            palette = viridis::magma(10))
#' }
#'
#' @export


plot_and_average <- function(seurat_obj, markers, palette = magma(10),
                             min_cutoff = "q10", max_cutoff = "q90", assay = "RNA",
                             reduction = "umap") {
  
  # loading packages
  list.of.packages = c('ggplot2', "viridis", "viridisLite" , "Seurat",'dplyr', 
                       'cowplot', "patchwork")
  
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
  if(length(new.packages) > 0) install.packages(new.packages)
  
  invisible(lapply(list.of.packages, FUN = library, character.only = T))
  rm(list.of.packages, new.packages)
  
  # 1. Graficar expresión de los marcadores
  ncol <- case_when(
    length(markers) == 1 ~ 1,
    length(markers) > 2  ~ 3,
    TRUE                 ~ 2
  )
  
  
  p <- FeaturePlot(seurat_obj,
                   features = markers,
                   ncol = ncol,
                   order = TRUE,
                   reduction = reduction,
                   #repel = TRUE,
                   label = TRUE,
                   cols = palette,
                   min.cutoff = min_cutoff,
                   max.cutoff = max_cutoff)
  
  print(p)  # Mostrar el gráfico
  
  # 2. Calcular expresión promedio
  avg_exp <- AverageExpression(seurat_obj, features = markers, assay = assay)[[assay]]
  
  # 3. Convertir a data frame y añadir columna con la media
  avg_exp_df <- avg_exp %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    rowwise() %>%
    mutate(media = mean(c_across(where(is.numeric)))) %>%
    as.data.frame()
  
  return(avg_exp_df)
}