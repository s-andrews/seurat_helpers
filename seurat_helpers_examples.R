library(Seurat)
library(ggpubr)
load("test_data.Rda")
source("seurat_helpers_functions.R")


# Additional QC metrics
add_qc_metrics(data) -> data
head(data[[]])


# Pseudobulk
c(
  "AAACCTGAGAAACCAT-1",
  "AAACCTGTCCAGATCA-1",
  "AAAGATGAGGAGCGTT-1",
  "AAAGATGCATTAACCG-1",
  "AAAGCAATCCAGTATG-1",
  "AAATGCCAGTGAACGC-1", 
  "AAATGCCCAATGGAAT-1"
  ) -> cell_ids1

c("AAACCTGAGATAGCAT-1",
  "AAACCTGCAACGCACC-1",
  "AAACGGGAGCTACCTA-1",
  "AAACGGGCAACTGCTA-1",
  "AAAGCAACATTAGGCT-1",
  "AACACGTGTTCGTCTC-1"
  ) -> cell_ids2

tibble(
  gene = rownames(data),
  group1 = pseudobulk(data,cell_ids1),
  group2 = pseudobulk(data,cell_ids2)
) -> bulk_data

bulk_data %>%
  ggplot(aes(x=group1+1,y=group2+1)) +
  geom_point() +
  scale_x_log10() + 
  scale_y_log10()


# Complexity
calculate_complexity(data) -> data

plot_complexity(data, limit=0.1)

# QC plots
plot_combined_qc(data) -> combined_qc_plots

ggarrange(plotlist = combined_qc_plots, ncol = 1)

# QC metrics per cluster
plot_cluster_qc(data) -> plot_cluster_plots

ggarrange(plotlist = plot_cluster_plots, ncol=1)

FindVariableFeatures(data) -> data
ScaleData(data) -> data
RunPCA(data) -> data

# QC metrics on PCA
plot_reduction_qc(data) -> reduction_qc_plots
ggarrange(plotlist = reduction_qc_plots, nrow=4, ncol=2)

# QC metrics on PCA log
plot_reduction_qc(data, log_scale = TRUE) -> reduction_qc_plots_log
ggarrange(plotlist = reduction_qc_plots_log, nrow=4, ncol=2)

# Knee plot
knee_plot(data)

