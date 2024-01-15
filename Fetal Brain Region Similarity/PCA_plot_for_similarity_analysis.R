
library(Seurat)
library(dplyr)
library(ggplot2)

# Assuming you have two Seurat objects: seurat_obj_your_data and seurat_obj_fetal_brain
load("~/Desktop/wk4_map2_slc32a1_only_mito10_stress_filt_CC_regressed_w_pyscenicAUC_GW20_annotation_20231114.RData")
load("~/Desktop/GW14_all_regions_merged_w_celltype_20231115.RData")
load("~/Documents/Ruiqi_scRNAseq/RData_files/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_neurons_only_badhuri_060802023.RData")
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = GAD1 > 0 & GAD2 > 0)
GW14_neurons <- subset(x = GW14, subset = (celltype == "Neuron" | celltype == "Interneuron" ))
save(GW14_neurons, file = "~/Desktop/GW14_all_regions_merged_w_celltype_neurons_only_20231117.RData")
GW14_neurons_GADexp_neurons_interneurons <- subset(GW14_neurons, subset = GAD1 > 0 & GAD2 > 0)
wk4_fetal_brain_seurat_integrated_noGW14 <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW14", invert = TRUE)
fetal_brain_seurat_integrated <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "wk4", invert = TRUE)
fetal_brain_seurat_integrated_GW20_34_only <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW20_34")
fetal_brain_seurat_integrated_GW20_new_only <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW20")
wk4_fetal_brain_seurat_integrated_GW20_34_wk4_only <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW20_34" | dataset == "wk4")
wk4_fetal_brain_seurat_integrated_GW20_new_wk4_only <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW20" | dataset == "wk4")
wk4_fetal_brain_seurat_integrated_wk4_only <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "wk4")

load("~/Desktop/GW20_new_individual_merged_regions_w_celltype_neurons_only_20231117.RData")
GW20_new_indvidual_neurons_GADexp_neurons_interneurons <- subset(GW20, subset = GAD1 > 0 & GAD2 > 0)


# Average expression for each gene across all cells in a cluster
wk4_pseudobulk <- AverageExpression(wk4_map2_slc32a1_only, assay = "RNA", group.by = "seurat_clusters")
# Convert to matrix
wk4_pseudobulk <- as.matrix(wk4_pseudobulk$RNA) 
write.csv(wk4_pseudobulk, file="wk4_pseudobulk.csv")


# Average expression for each gene across all cells in a cluster
GW20_new_indiv_pseudobulk <- AverageExpression(fetal_brain_seurat_integrated_GW20_new_only, assay = "RNA", group.by = "structure")
GW20_new_indiv_pseudobulk <- as.matrix(GW20_new_indiv_pseudobulk$RNA) 

GW20_34_indiv_pseudobulk <- AverageExpression(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, assay = "RNA", group.by = "structure")
GW20_34_indiv_pseudobulk <- as.matrix(GW20_34_indiv_pseudobulk$RNA) 

wk4_fetal_brain_seurat_integrated_pseudobulk <- AverageExpression(wk4_fetal_brain_seurat_integrated_GW20_new_wk4_only, assay = "RNA", group.by = "structure")
wk4_fetal_brain_integrated_pseudobulk_pseudobulk <- as.matrix(wk4_fetal_brain_seurat_integrated_pseudobulk$RNA) 

wk4_rows <- rownames(wk4_pseudobulk)
GW20_34_rows <- rownames(GW20_34_indiv_pseudobulk)
GW20_rows <- rownames(GW20_new_indiv_pseudobulk)
wk4_fetal_brain_integrated_pseudobulk_pseudobulk_rows <- rownames(wk4_fetal_brain_integrated_pseudobulk_pseudobulk)


common_genes <- intersect(wk4_fetal_brain_integrated_pseudobulk_pseudobulk_rows, GW20_34_rows)
common_genes_all <- intersect(common_genes, GW20_rows)
dim(common_genes_all)
var_genes <- wk4_fetal_brain_seurat_integrated@assays$integrated@var.features
wk4_fetal_brain_seurat_integrated_GW20_new_wk4_only <- FindVariableFeatures(wk4_fetal_brain_seurat_integrated_GW20_new_wk4_only, nfeatures = 2000)
var_genes <- wk4_fetal_brain_seurat_integrated_GW20_new_wk4_only@assays$RNA@var.features

DefaultAssay(wk4_fetal_brain_seurat_integrated) <- "RNA"

# Subset datasets 
wk4_sub <- wk4_pseudobulk[var_genes,]
GW20_34_sub <- GW20_34_indiv_pseudobulk[var_genes,]
GW20_sub <- GW20_new_indiv_pseudobulk[var_genes,]
wk4_fetal_brain_integrated_pseudobulk_sub <- wk4_fetal_brain_integrated_pseudobulk_pseudobulk[var_genes,]

# Combine datasets
combined_data <- cbind(wk4_fetal_brain_integrated_pseudobulk_sub, GW20_34_sub)


# Normalize 
combined_matrix_norm <- apply(wk4_fetal_brain_integrated_pseudobulk_sub, 2, function(x) (x - mean(x)) / sd(x))

# Scale 
combined_matrix_scaled <- scale(combined_matrix_norm)

# Run PCA
pca <- prcomp(combined_matrix_scaled, scale = FALSE) 

# Extract PC vectors 
#pca_data <- pca$x

pca_rotation <- pca$rotation

# Convert to data frame
pca_rotation_df <- as.data.frame(pca_rotation)

library(ggfortify)
library(ggrepel)
ggbiplot(pca_rotation_df, obs.scale = 1, var.scale = 1, 
         groups = colnames(combined_matrix_scaled), ellipse = TRUE) +
  ggrepel::geom_label_repel(
    aes(x=PC1, y=PC2, label=rownames(pca_rotation_df))
  )


 # Make scatter plot
#plot(pca_loadings[,1], pca_loadings[,2], 
    # col=c(rep("blue", nrow(GW20_34_pseudobulk)), rep("blue", nrow(wk4_pseudobulk))),
    # xlab = "PC1", ylab = "PC2")



#method 2
# For Seurat object 1
wk4_pseudobulk <- AverageExpression(wk4_map2_slc32a1_only, assay = "RNA", group.by = "seurat_clusters")
wk4_pseudobulk_data <- as.data.frame(wk4_pseudobulk$RNA)


# For Seurat object 2
GW20_34_pseudobulk <- AverageExpression(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, assay = "RNA", group.by = "structure")
GW20_34_pseudobulk_data <- as.data.frame(GW20_34_pseudobulk$RNA)

GW14_pseudobulk <- AverageExpression(GW14_neurons_GADexp_neurons_interneurons, assay = "RNA", group.by = "structure")
GW14_pseudobulk_data <- as.data.frame(GW14_pseudobulk$RNA)

GW20_new_pseudobulk <- AverageExpression(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, assay = "RNA", group.by = "structure")
GW20_new_pseudobulk_data <- as.data.frame(GW20_new_pseudobulk$RNA)

#merge
wk4_pseudobulk_data$dataset <- "Wk4"
GW20_34_pseudobulk_data$dataset <- "GW20_34"
GW20_new_pseudobulk_data$dataset <- "GW20_new"
GW14_pseudobulk_data$dataset <- "GW14"

# Assuming df1 and df2 are your data frames
wk4_pseudobulk_data$genes <- rownames(wk4_pseudobulk_data)
GW20_34_pseudobulk_data$genes <- rownames(GW20_34_pseudobulk_data)
GW20_new_pseudobulk_data$genes <- rownames(GW20_new_pseudobulk_data)
GW14_pseudobulk_data$genes <- rownames(GW14_pseudobulk_data)


library(dplyr)

combined_df <- wk4_pseudobulk_data %>%
  full_join(GW20_34_pseudobulk_data, by = "dataset") %>%
  full_join(GW20_new_pseudobulk_data, by = "dataset") %>%
  full_join(GW14_pseudobulk_data, by = "dataset") 

combined_df <- merge(wk4_pseudobulk_data, GW20_34_pseudobulk_data, by = "genes")

# Ensure that the data is numeric for PCA
combined_pseudobulk_numeric <- combined_df[, sapply(combined_df, is.numeric)]

# Perform PCA
pca_result <- prcomp(combined_pseudobulk_numeric, center = TRUE, scale. = TRUE)

library(ggplot2)

pca_rotation <- pca_result$rotation

# Convert to data frame
pca_rotation_df <- as.data.frame(pca_rotation)

library(ggfortify)
library(ggrepel)
ggbiplot(pca_rotation_df, obs.scale = 1, var.scale = 1, 
         groups = colnames(combined_matrix_scaled), ellipse = TRUE) +
  ggrepel::geom_label_repel(
    aes(x=PC1, y=PC2, label=rownames(pca_rotation_df))
  )





#############
#chnage headings of seurat objects prior to merge to include age for PCA plot
# Assuming your Seurat object is named 'seurat_obj'
# And the column you want to modify is named 'your_column_name'

# Access and modify the column
metadata <- GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons@meta.data
structure <- metadata$structure
test <- paste0("GW20_", structure)
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- AddMetaData(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, metadata = test, col.name = "structure")
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- AddMetaData(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, metadata = "GW20_34", col.name = "dataset")
head(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons)

metadata <- GW20_new_indvidual_neurons_GADexp_neurons_interneurons@meta.data
structure <- metadata$structure
test <- paste0("GW20_", structure)
GW20_new_indvidual_neurons_GADexp_neurons_interneurons <- AddMetaData(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, metadata = test, col.name = "structure")
GW20_new_indvidual_neurons_GADexp_neurons_interneurons <- AddMetaData(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, metadata = "GW20", col.name = "dataset")
head(GW20_new_indvidual_neurons_GADexp_neurons_interneurons)

metadata <- wk4_map2_slc32a1_only@meta.data
structure <- metadata$seurat_clusters
test <- paste0("Wk4", structure)
wk4_map2_slc32a1_only <- AddMetaData(wk4_map2_slc32a1_only, metadata = structure, col.name = "structure")
wk4_map2_slc32a1_only <- AddMetaData(wk4_map2_slc32a1_only, metadata = "wk4", col.name = "dataset")
head(wk4_map2_slc32a1_only)

metadata <- GW14_neurons_GADexp_neurons_interneurons@meta.data
structure <- metadata$structure
test <- paste0("GW14", structure)
GW14_neurons_GADexp_neurons_interneurons <- AddMetaData(GW14_neurons_GADexp_neurons_interneurons, metadata = test, col.name = "structure")
GW14_neurons_GADexp_neurons_interneurons <- AddMetaData(GW14_neurons_GADexp_neurons_interneurons, metadata = "GW14", col.name = "dataset")
head(GW14_neurons_GADexp_neurons_interneurons)


objects <- ls(pat = "GW14")
object <- objects[-c(1,3,8)]
cell_ids <- c("CGE", "hypo", "LGE", "MGE", "motor", "occipital", "somato", "striatum", "thalamus")

#merge
combined <- merge(wk4_map2_slc32a1_only, y = c(GW14_neurons_GADexp_neurons_interneurons, GW20_new_indvidual_neurons_GADexp_neurons_interneurons, GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons), add.cell.ids = c("wk4", "GW14", "GW20_34", "GW20"),  project = "brain_similarity")

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

save(combined, file = "~/Desktop/wk4_merged_w_GW14_GW20_GW20_34_fetal_brain_regions_20231117.RData")

#integration of striatum and hypo datasets

# split the merged  dataset into a list  seurat objects 
ifnb.list <- SplitObject(combined, split.by = "dataset")


# normalize and identify variable features for each dataset independently
ifnb.list <- c(GW20_caudate, GW20_hypo, GW20_NA, GW20_putamen)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list, dims = 1:30)

#perform integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
wk4_fetal_brain_seurat_integrated <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(wk4_fetal_brain_seurat_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
wk4_fetal_brain_seurat_integrated <- ScaleData(wk4_fetal_brain_seurat_integrated, verbose = FALSE)
wk4_fetal_brain_seurat_integrated <- FindVariableFeatures(wk4_fetal_brain_seurat_integrated, nfeatures = 2000)
wk4_fetal_brain_seurat_integrated <- RunPCA(wk4_fetal_brain_seurat_integrated, npcs = 30, dims = 1:30)
wk4_fetal_brain_seurat_integrated <- RunUMAP(wk4_fetal_brain_seurat_integrated, reduction = "pca", dims = 1:30)
wk4_fetal_brain_seurat_integrated <- RunUMAP(wk4_fetal_brain_seurat_integrated, dims = 1:30)
wk4_fetal_brain_seurat_integrated <- FindNeighbors(wk4_fetal_brain_seurat_integrated)
wk4_fetal_brain_seurat_integrated <- FindClusters(wk4_fetal_brain_seurat_integrated, resolution = 0.6)

# Visualization
DimPlot(wk4_fetal_brain_seurat_integrated_noGW14, reduction = "pca", group.by = "structure")
DimPlot(wk4_fetal_brain_seurat_integrated_noGW14, reduction = "umap", group.by = "structure", pt.size = .001)
p2 <- DimPlot(gcdata_seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


DimPlot(object = wk4_fetal_brain_seurat_integrated_noGW14, dims = c(2, 5), reduction = "pca", group.by = "structure")

save(wk4_fetal_brain_seurat_integrated, file = "~/Desktop/wk4_fetal_brain_seurat_integrated_20231117.RData")


# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(GW20_hypo_stri_regions_combined, reduction = "umap", split.by = "stim")


#from chatGPT for pearson/spearman correlation 
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)



# Identify clusters or cell types - assuming this is already done
# seurat_object1 <- FindClusters(seurat_object1, resolution = 0.5)
# seurat_object2 <- FindClusters(seurat_object2, resolution = 0.5)

wk4_from_integration <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "wk4")
wk4_fetal_brain_seurat_integrated_GW20_34_only <- subset(wk4_fetal_brain_seurat_integrated_noGW14, subset = dataset == "GW20_34")

# Pseudobulk approach - aggregate data by cluster
pseudo_bulk1 <- AverageExpression(wk4_from_integration, group.by = 'structure', return.seurat = TRUE)
pseudo_bulk2 <- AverageExpression(wk4_fetal_brain_seurat_integrated_GW20_34_only, group.by = "structure", return.seurat = TRUE)


# Extract averaged data
avg_expr1 <- pseudo_bulk1@assays$RNA@data
avg_expr2 <- pseudo_bulk2@assays$RNA@data

# Ensure same set of genes in both matrices for correlation
common_genes <- intersect(rownames(avg_expr1), rownames(avg_expr2))
use_genes <- intersect(common_genes, var_genes)
avg_expr1 <- avg_expr1[var_genes, ]
avg_expr2 <- avg_expr2[var_genes, ]

# Calculate Pearson correlation
pearson_corr <- cor(avg_expr1, avg_expr2, method = "pearson")
write.csv(pearson_corr, file = "~/Documents/Ruiqi_scRNAseq/wk4_GW20_perason_correlation_values_for_heatmap_20231122.csv")

# Calculate Spearman correlation
spearman_corr <- cor(avg_expr1, avg_expr2, method = "spearman")

# Heatmap Visualization for Pearson Correlation
pheatmap(pearson_corr, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Pearson Correlation Heatmap")

# Heatmap Visualization for Spearman Correlation
pheatmap(spearman_corr, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Spearman Correlation Heatmap")



