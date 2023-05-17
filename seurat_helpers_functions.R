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
  # For this we'll remove all mitochondria, ribosomal and malat data
  genes_to_remove = str_detect(tolower(rownames(data)),"^mt-") | str_detect(tolower(rownames(data)),"^rp[ls]") | tolower(rownames(data))=="malat1"
  
  data@assays$RNA@counts[!genes_to_remove,] -> data.nomalat
  
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

largest_gene_table <- function (data) {
  data[[]] %>%
    group_by(Largest_Gene) %>%
    count(name = "number_of_cells") %>%
    arrange(desc(number_of_cells)) %>%
    ungroup() %>%
    return()
}

largest_gene_per_cluster_table <- function (data) {
  filtered_data[[]] %>%
    group_by(Largest_Gene,seurat_clusters) %>%
    count(name = "number_of_cells") %>%
    arrange(desc(number_of_cells)) %>%
    group_by(seurat_clusters) %>%
    slice(1:5) %>%
    ungroup() %>%
    arrange(seurat_clusters,desc(number_of_cells)) %>%
    return()
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


#' Plot combinations of QC metrics
#'
#' @param data a Seurat object
#'
#' @return A list of 3 plots
#' @export
#'
#' @examples
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
    slice(1:10) %>%
    pull(Largest_Gene) -> largest_genes_to_plot
  
  data[[]] %>%
    filter(Largest_Gene %in% largest_genes_to_plot) %>%
    mutate(Largest_Gene=factor(Largest_Gene, levels=largest_genes_to_plot)) %>%
    arrange(Largest_Gene) %>%
    ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=Largest_Gene)) +
    geom_point(size=1) +
    scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1"))) + 
    guides(colour = guide_legend(override.aes = list(size=4))) -> plot2
  
  data[[]] %>%
    filter(Largest_Gene %in% largest_genes_to_plot) %>%
    mutate(Largest_Gene=factor(Largest_Gene, levels=largest_genes_to_plot)) %>%
    arrange(Largest_Gene) %>%
    ggplot(aes(x=complexity, y=percent_Largest_Gene, colour=Largest_Gene)) +
    geom_point()+
    scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1"))) + 
    guides(colour = guide_legend(override.aes = list(size=4))) -> plot3
  
  
  return(list(plot1,plot2,plot3))
  
}


#' Plot all QC metrics for all clusters
#'
#' @param data a Seurat object
#'
#' @return a list of two plots
#' @export
#'
#' @examples
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
    geom_boxplot(show.legend = FALSE) +
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



#' Plot all QC metrics on a dimension reduction plot
#'
#' @param data the seurat object
#' @param reduction which reduction to use (default pca)
#' @param dims which dimensions to use
#' @param desc reverse the order of the values
#' @param log_scale plot on a log scale
#'
#' @return a list of plots
#' @export
#'
#' @examples
plot_reduction_qc <- function(data,reduction="pca", dims=c(1,2), desc=FALSE, log_scale=FALSE) {
  
  
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
    if (log_scale) {
      plot_data %>%
        ggplot(aes_string(x=embedding_names[1], y=embedding_names[2])) +
        geom_point(aes(colour=log10(!!sym(metric)))) +
        scale_colour_gradientn(colours = c("magenta3","grey","green3")) -> p
      
    }else {
      plot_data %>%
        ggplot(aes_string(x=embedding_names[1], y=embedding_names[2], colour=metric)) +
        geom_point() +
        scale_colour_gradientn(colours = c("magenta3","grey","green3")) -> p
      
    }
    
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

#' Plot a knee plot of your data
#'
#' @param data The seurat object to plot
#'
#' @return a knee plot of the data
#' @export
#'
#' @examples
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

#' Plot integration anchors
#'
#' @param data_split The list of split datasets to integrate
#' @param data The original data containing the tsne embedding
#' @param anchors The set of defined anchors from FindIntegrationAnchors
#' @param data1 The index of the first dataset to plot
#' @param data2 The index of the second dataset to plot
#' @param min_score The minimum anchor score to plot
#'
#' @return A ggplot of the anchors
#' @export
#'
#' @examples
plot_anchors <- function(data_split, data, anchors, data1,data2,min_score=0.5) {
  
  colnames(data_split[[data1]]) -> c1
  colnames(data_split[[data2]]) -> c2
  
  tibble(
    cell_id=c1,
    cell=1:length(c1)
  ) -> c1
  
  tibble(
    cell_id=c2,
    cell=1:length(c2)
  ) -> c2
  
  anchors@anchors %>%
    filter(dataset1==data1 & dataset2==data2)  %>%
    filter(score>min_score) %>%
    mutate(pair=1:n()) -> anchor_cells
  
  anchor_cells %>%
    select(pair,cell1) %>%
    rename(cell=cell1) %>%
    left_join(c1) %>%
    select(-cell) -> c1
  
  anchor_cells %>%
    select(pair,cell2) %>%
    rename(cell=cell2) %>%
    left_join(c2) %>%
    select(-cell) -> c2
  
  bind_rows(c1,c2) %>%
    arrange(pair) -> paired_cells
  
  rm(c1, c2)
  
  data@reductions$tsne@cell.embeddings %>%
    as_tibble(rownames="cell_id") %>%
    right_join(paired_cells) %>%
    select(-cell_id) %>%
    arrange(pair) %>%
    mutate(from_to=rep(c("from","to"),n()/2)) %>%
    pivot_wider(
      names_from=from_to,
      values_from=c(tSNE_1,tSNE_2)
    ) -> paired_tsne
  
  colnames(anchors@object.list[[1]])
  
  names(anchors@object.list)[1]
  
  data@reductions$tsne@cell.embeddings %>%
    as_tibble(rownames="cell") %>%
    mutate(dataset=case_when(
      cell %in% colnames(anchors@object.list[[data1]]) ~   names(anchors@object.list)[data1],
      cell %in% colnames(anchors@object.list[[data2]]) ~   names(anchors@object.list)[data2]
    )) %>%
    filter(!is.na(dataset)) %>%
    ggplot(aes(x=tSNE_1, y=tSNE_2, colour=dataset)) +
    geom_point() +
    scale_colour_brewer(palette = "Set1") +
    geom_segment(data=paired_tsne, aes(x=tSNE_1_from, y=tSNE_2_from, xend=tSNE_1_to, yend=tSNE_2_to), colour="black") +
    xlab(data1) + ylab(data2) -> p
  
  return(p)
  
}













