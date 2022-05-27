library(Seurat)
load("test_data.Rda")

add_qc_metrics(data) -> data
head(data[[]])