
#' Adds some extra QC metrics to a seurat data object
#'
#' @param data A seurat data object
#'
#' @return A seurat data object with additional QC metrics
#' @export
#'
#' @examples
add_qc_metrics <- function(data) {
  
  # Mitochrondrial sequences
  message("Adding percent mitochondrial")
  PercentageFeatureSet(data,pattern="^MT-") -> data$percent_MT
  
  # Ribosomal sequences
  message("Adding percent ribosomal")
  PercentageFeatureSet(data,pattern="^RP[LS]") -> data$percent_Ribosomal
  
  # Malat1
  message("Adding percent Malat1")
  PercentageFeatureSet(data,pattern="MALAT1") -> data$percent_Malat
  
  # Largest Gene
  message("Adding largest gene information")
  
  message("Making data subset")
  data@assays$RNA@counts[rownames(data) != "MALAT1",] -> data.nomalat
  
  message("Getting largest count")
  apply(
    data.nomalat,
    2,
    max
  ) -> largest_count
  
  message("Finding largest index")
  apply(
    data.nomalat,
    2,
    which.max
  ) -> largest_index
  
  message("Getting largest name")
  rownames(data)[largest_index] -> data$Largest_Gene
  
  message("Calculating percent largest")
  100 * largest_count / data$nCount_RNA -> data$percent_Largest_Gene
  
  rm(data.nomalat)
  
  return(data)
  
}

calculate_complexity <- function(data) {
  log10(data$nFeature_RNA) / log10(data$nCount_RNA)  -> complexity
  
  lm(log10(data$nFeature_RNA)~log10(data$nCount_RNA)) -> complexity.lm
  
  data$complexity = log10(data$nFeature_RNA) - ((log10(data$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
  
  return(data)
}




#' Get Pseudobulk counts for a collection of cells
#'
#' @param data A seurat data object
#' @param cell_ids A vector of cell ids to group together
#'
#' @return A vector of 
#' @export
#'
#' @examples
pseudobulk <- function(data, cell_ids) {
  
  return(rowSums(data@assays$RNA@counts[,cell_ids]))
  
}




