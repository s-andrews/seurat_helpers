library(Seurat)
library(tidyverse)


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
  if (any(startsWith(rownames(data),"mt-"))) {
    PercentageFeatureSet(data,pattern="^mt-") -> data$percent_MT
  }else {
    PercentageFeatureSet(data,pattern="^MT-") -> data$percent_MT
  }
  
  # Ribosomal sequences
  message("Adding percent ribosomal")
  if (any(startsWith(rownames(data),"Rpl"))) {
    PercentageFeatureSet(data,pattern="^Rp[ls]") -> data$percent_Ribosomal
  }else {
    PercentageFeatureSet(data,pattern="^RP[LS]") -> data$percent_Ribosomal
  }

  # Malat1
  message("Adding percent Malat1")
  if ("MALAT1" %in% rownames(data)) {
    PercentageFeatureSet(data,pattern="MALAT1") -> data$percent_Malat
  }else {
    PercentageFeatureSet(data,pattern="Malat1") -> data$percent_Malat
  }
  
  # Largest Gene
  message("Adding largest gene information")
  
  message("Making data subset")
  data@assays$RNA@counts[!rownames(data) %in% c("MALAT1","Malat1"),] -> data.nomalat
  
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


#' Plot out the complexity values for a cell population
#'
#' @param data The seurat object holding data for which calculate_complexity has been run
#' @param limit The limit for the colour scale on the plot
#'
#' @return The ggplot graph of the plot
#' @export
#'
#' @examples
plot_complexity <- function(data, limit=0.1) {
  as_tibble(data[[]]) %>%
    mutate(complexity=replace(complexity,complexity < -limit, -limit)) %>%
    mutate(complexity=replace(complexity,complexity > limit, limit)) %>%
    ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity)) +
    geom_point(size=1) +
    scale_colour_gradient2(low="blue2",mid="grey",high="red2") %>%
    return()
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


plot_combined_qc <- function (data) {
  
  # Counts vs features scatterplot
  data[[]] %>%
    arrange(percent_MT) %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, colour=percent_MT)) +
    geom_point(size=1) +
    scale_x_log10() +
    scale_y_log10() + 
    scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) -> plot1
  
  
  # Complexity vs largest gene
  data[[]] %>%
    group_by(Largest_Gene) %>%
    count() %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    slice(1:9) %>%
    pull(Largest_Gene) -> largest_genes_to_plot
  
  data[[]] %>%
    filter(Largest_Gene %in% largest_genes_to_plot) %>%
    mutate(Largest_Gene=factor(Largest_Gene, levels=largest_genes_to_plot)) %>%
    arrange(Largest_Gene) %>%
    ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=Largest_Gene)) +
    geom_point(size=1) +
    scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1"))) -> plot2
  
  data[[]] %>%
    filter(Largest_Gene %in% largest_genes_to_plot) %>%
    mutate(Largest_Gene=factor(Largest_Gene, levels=largest_genes_to_plot)) %>%
    arrange(Largest_Gene) %>%
    ggplot(aes(x=complexity, y=percent_Largest_Gene, colour=Largest_Gene)) +
    geom_point()+
    scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1"))) -> plot3
  
  
  return(list(plot1,plot2,plot3))
  
}


plot_cluster_qc <- function(data) {
  
  data[[]] %>%
    as_tibble() %>%
    filter(!is.na(seurat_clusters)) %>%
    select(seurat_clusters,where(is.numeric)) %>%
    pivot_longer(
      cols=-seurat_clusters,
      names_to="QC_Metric",
      values_to="Value"
    ) -> plot_data
  
  
  plot_data %>%
    ggplot(aes(x=seurat_clusters, y=Value, fill=QC_Metric)) +
    geom_violin(show.legend = FALSE) +
    facet_wrap(vars(QC_Metric), scale="free_y") -> p
  
  plot_data %>%
    group_by(QC_Metric) %>%
    mutate(Z_Value=(Value-mean(Value))/sd(Value)) %>%
    mutate(Z_Value=replace(Z_Value, Z_Value>5,5)) %>%
    mutate(Z_Value=replace(Z_Value, Z_Value< -5,-5)) %>%
    ggplot(aes(x=QC_Metric, y=Z_Value, fill=QC_Metric)) +
    geom_violin(show.legend = FALSE) +
    geom_hline(yintercept = 0, size=1) +
    coord_flip() +
    facet_wrap(vars(seurat_clusters)) -> p2
  
  return(list(p,p2))      
  
}


plot_reduction_qc <- function(data,reduction="pca", dims=c(1,2), desc=FALSE) {
  
  
  # Extract the relevant embeddings
  Embeddings(data, reduction=reduction)[,dims] %>%
    as_tibble(rownames="cell_id") -> plot_data
  
  colnames(plot_data)[2:3] -> embedding_names
  
  data[[]] %>%
    as_tibble(rownames="cell_id") %>%
    select(cell_id,where(is.numeric)) %>%
    right_join(plot_data) -> plot_data
  
  
  data[[]] %>%
    as_tibble(rownames="cell_id") %>%
    select(where(is.numeric)) %>%
    colnames() -> plot_metrics
  
  
  lapply(plot_metrics, function(metric){
    if(desc) {
      plot_data %>%
        arrange(desc(!!sym(metric))) -> plot_data
    }else{
      plot_data %>%
        arrange(!!sym(metric)) -> plot_data
    }
    plot_data %>%
      ggplot(aes_string(x=embedding_names[1], y=embedding_names[2], colour=metric)) +
      geom_point() +
      scale_colour_gradientn(colours = c("magenta3","grey","green3")) -> p
    
    return(p)
    
  }) -> return_list
  

  return(return_list)      
  
}


plot_largest_gene_dimensions <- function(data, reduction="pca", dims=c(1,2)) {
  
  # Extract the relevant embeddings
  Embeddings(data, reduction=reduction)[,dims] %>%
    as_tibble(rownames="cell_id") -> plot_data
  
  Embeddings(data, reduction=reduction)[,dims] %>%
    colnames() -> embedding_names
  
  data[[]] %>%
    as_tibble(rownames="cell_id") %>%
    select(cell_id,seurat_clusters,Largest_Gene) %>%
    right_join(plot_data) %>%
    filter(!is.na(seurat_clusters)) -> plot_data

  plot_data %>%
    group_by(seurat_clusters,Largest_Gene) %>%
    count() %>%
    arrange(desc(n)) %>%
    group_by(seurat_clusters) %>%
    slice(1) %>%
    ungroup() %>%
    select(-n) %>%
    rename(Largest_Gene_Per_Cluster = Largest_Gene) %>%
    right_join(plot_data) -> plot_data
  
  
  plot_data %>%
    ggplot(aes_string(x=embedding_names[1], y=embedding_names[2], colour="Largest_Gene_Per_Cluster")) +
    geom_point() +
    scale_colour_brewer(palette = "Set1") + 
    guides(colour = guide_legend(override.aes = list(size=4))) -> p

  
  plot_data %>%
    ggplot(aes_string(x=embedding_names[1], y=embedding_names[2], colour="Largest_Gene_Per_Cluster")) +
    geom_point() +
    scale_colour_brewer(palette = "Set1") + 
    guides(colour = guide_legend(override.aes = list(size=4))) +
    facet_wrap(vars(seurat_clusters)) -> p2
  
  
    
  return(list(p,p2))
  
}

knee_plot <- function(data) {
  data[[]] %>%
    as_tibble() %>%
    select(nCount_RNA) %>%
    arrange(desc(nCount_RNA)) %>%
    mutate(index=1:n()) %>%
    ggplot(aes(x=index, y=nCount_RNA)) +
    geom_line() +
    scale_x_log10() +
    scale_y_log10() +
    xlab("Number of cells") +
    ylab("Number of UMIs") -> p
  
  return(p)
}













