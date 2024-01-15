install.packages('devtools')
install.packages("cli")
library(devtools)

devtools::install_github('immunogenomics/presto')

devtools::install_github('quadbiolab/voxhunt')

#for organoid data using seurat
# Load VoxHunt
library(voxhunt)
library(patchwork)
library(tidyverse)

#load human brainspan data
data('brainspan')
names(brainspan)

voxhunt::plot_annotation('8 pcw')

# Point VoxHunt to ABA expression data
load_aba_data('~/Downloads/voxhunt_data/')

# Find 300 most variable genes from the E13.5 mouse brain
genes_use <- variable_genes('E13', 300)$gene

# Calculate the similarity map of a seurat object to the E13.5 mouse brain
vox_map <- voxel_map(wk4_map2_slc32a1_only, genes_use=genes_use)

# Plot the result
plot_map(vox_map)


#using brainspan data for bulk RNAseq
regional_markers <- structure_markers('E15') %>%
  group_by(group) %>%
  top_n(10, auc) %>% 
  {unique(.$gene)}

#


ref_map <- brainspan_map(
  wk4_map2_slc32a1_only,
  stages = 10:24,
  genes_use = regional_markers,
  group_name = "seurat_clusters",
  method = 'pearson',
  pseudobulk_groups = T
)

print(ref_map)

plot_map(ref_map)
plot_structure_similarity(ref_map, annotation_level = 'structure', scale = F)


ref_map <- reference_map(
  wk4_map2_slc32a1_only,
  fetal_brain_integrated,
  slot = "data",
  assay = "RNA",
  group_name = 'seurat_clusters',
  reduction = "umap",
  method = "pearson",
  genes_use = g,
  allow_neg = FALSE,
  pseudobulk_groups = TRUE
)

