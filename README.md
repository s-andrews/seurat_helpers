# Seurat Helpers

This is a collection of additional R functions which can be used as part of an analysis using the Seurat package.

## Installation
At the moment this is just a bare R file of functions.  Once it's been developed more I will probably put it into a proper package structure.  For now you need to clone this git repository to your system.

```
git clone https://github.com/s-andrews/seurat_helpers.git
```

To use this in your script you'll then need to run the main file of functions:

```
source("seurat_helpers/seurat_helpers_functions.R")
```

You'll need to modify this path to be the absolute path to your git repository, or set your working directory appropriately so R can find the file.

## Usage

### QC

```add_qc_metrics```

This function puts some additional QC metrics into the main metadata table on a seurat object.  The current metrics which are added are:

* Percent Mitochrondrial Reads
* Percent Ribosomal Protein Reads
* Percent Malat1 reads
* The name of the most highly expressed (non-Malat1) gene
* The percentage of reads coming from the highest expressed gene


```calculate_complexity```

This function calculates a value to indicate the relative complexity of the cells in the library.  This is a measure which relates the number of reads to the number of genes detected.  Libraries with lower complexity might be dominated by the activity of a single, or small group of, genes.  The new value is added to the metadata table.

```plot_complexity```

This is a plot which shows the number of reads vs the number of genes detected, and is coloured by the complexity value created by ```calcuate_complexity```

```plot_combined_qc```
Generates plots of combinations of QC parameters which are often helpful to view in combination.

1. Counts vs Features coloured by percent MT

2. Counts vs Features coloured by largest gene

3. Complexity vs Percent Largest Gene


```plot_cluster_qc```
Plots all QC metrics split by cluster as a violin plot

```plot_reduction_qc```
Plots all QC metrics superimposed on a dimension reduction plot.

```knee_plot```
Plots a knee plot of your data

```plot_integration_anchors```
Plots the position of integration anchors on two datasets being integrated


### Pseudobulk
```pseudobulk```

One option for viewing data across a large number of cells is to create a pseudobulk sample, where we sum the counts from a subset of cells to produce a quantitation with much higher levels of observation.  You'd normally do this from a collection of cells which fell into the same cluster in an analysis to look at the overall expression of the cluster.
