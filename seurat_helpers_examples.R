library(Seurat)
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

plot_complexity(data)



