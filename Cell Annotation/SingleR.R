
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")

library(SingleR)
library(Seurat)
load("~/Desktop/wk4_fetal_brain_seurat_integrated_gad_exp_only_20231121.RData")

# Convert Seurat object to SingleCellExperiment
wk4_sce <- as.SingleCellExperiment(wk4_map2_slc32a1_only, assay = "RNA")
assayNames(wk4_sce)[assayNames(wk4_sce) == "existingAssayName"] <- "logcounts"
# Create a new logcounts assay - assuming 'counts' is your raw count data
logcounts(wk4_sce) <- log1p(counts(wk4_sce))

library(SingleCellExperiment)

# Assuming you have a data frame `expr_matrix` with rows as genes and columns as samples
# and a vector `cell_types` with the cell type annotation for each column in `expr_matrix`

fetal_brain_seurat_integrated <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "wk4", invert = TRUE)
fetal_brain_seurat_integrated_noGW14 <- subset(fetal_brain_seurat_integrated, subset = dataset == "GW14", invert = TRUE)

expr_matrix <- fetal_brain_seurat_integrated@assays$RNA@counts


fetal_sce <- as.SingleCellExperiment(fetal_brain_seurat_integrated_noGW14, assay = "RNA")
fetal_sce <- logNormCounts(fetal_sce)
Idents(fetal_brain_seurat_integrated_noGW14) <- "structure"
labels <- Idents(fetal_brain_seurat_integrated_noGW14)
#use aggregate function to create pseudobulk. takes into account within label heterogeneity
aggr <- aggregateReference(fetal_sce, labels, power=0)


# Retrieve expression data
data <- GetAssayData(fetal_brain_seurat_integrated_noGW14, assay = "RNA", slot = "counts")  # Use "data" slot if you need normalized data

# Retrieve cluster information
Idents(fetal_brain_seurat_integrated_noGW14) <- "structure"
clusters <- Idents(fetal_brain_seurat_integrated_noGW14)

# Initialize an empty matrix for pseudobulk data
pseudobulk_data <- matrix(nrow = nrow(data), ncol = length(unique(clusters)))
rownames(pseudobulk_data) <- rownames(data)
colnames(pseudobulk_data) <- paste0(unique(clusters))

# Aggregate data by cluster
for (cluster in unique(clusters)) {
  cells_in_cluster <- WhichCells(fetal_brain_seurat_integrated_noGW14, idents = cluster)
  pseudobulk_data[, cluster] <- rowMeans(data[, cells_in_cluster])
}

# 'pseudobulk_data' now contains your pseudobulk expression data
head(pseudobulk_data)


# Create a SingleCellExperiment object
ref_sce <- SingleCellExperiment(assays = list(counts = pseudobulk_data))
#ref_sce <- SingleCellExperiment(assays = list(counts = as.matrix(pseudobulk_data)))
structure <- fetal_brain_seurat_integrated_noGW14$structure
metadata(ref_sce)$structure <- structure
# Rename an existing assay - replace 'existingAssayName' with the actual name of the assay
assayNames(ref_sce)[assayNames(ref_sce) == "existingAssayName"] <- "logcounts"
# Create a new logcounts assay - assuming 'counts' is your raw count data
library(scuttle)
ref_sce <- logNormCounts(ref_sce) 
#logcounts(ref_sce) <- log1p(counts(ref_sce))



# Assuming `sce` is your experimental single-cell data in SingleCellExperiment format
annotations <- SingleR(test = wk4_sce, ref = aggr, labels = colnames(pseudobulk_data))
annotations <- SingleR(test = wk4_sce, ref = aggr, labels = colnames(pseudobulk_data), clusters = wk4_map2_slc32a1_only$seurat_clusters)



# Load a built-in reference dataset
#library(AnnotationHub)
#ah <- AnnotationHub()
#ref <- query(ah, "SingleR")
#ref_dataset <- ref[['AHXXXXXX']]  # Replace 'AHXXXXXX' with the specific reference dataset ID

# Run SingleR
annotations <- SingleR(test = sce, ref = ref_dataset, clusters = wk4_map2_slc32a1_only$seurat_clusters)

# Add annotations to Seurat object
wk4_map2_slc32a1_only$SingleR_labels <- annotations$labels

DimPlot(wk4_map2_slc32a1_only, reduction = "umap", group.by = "SingleR_labels")

redblue<-grDevices::colorRampPalette(c('blue', 'white', 'red'))(100)
plotScoreHeatmap(annotations, color = redblue)

pheatmap(annotations$scores, cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 4, col = color, clustering_distance_rows = "correlation", clustering_method = "average")

classifySingleR(test, trained, quantile=0.8, fine.tune=TRUE, tune.thresh=0.05, 
                 sd.thresh=NULL, prune=TRUE, assay.type="logcounts", check.missing=TRUE, 
                 num.threads=bpnworkers(BPPARAM), BPPARAM=SerialParam() )

# Normalize 
combined_matrix_norm <- apply(annotations$scores, 2, function(x) (x - mean(x)) / sd(x))

# Scale 
combined_matrix_scaled <- scale(annotations$scores)

annotations$scores <- combined_matrix_scaled

#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#apply Min-Max normalization to first four columns in iris dataset
new_pseudo1 <- min_max_norm(annotations$scores)
annotations$scores <- new_pseudo1

om.all$cell_type_by_cluster <- factor(
  om.all$seurat_clusters,
  levels = rownames(cluster_annots)[cluster_anno_keep],
  labels = annots$pruned.labels[cluster_anno_keep])

seurat.obj[["SingleR.cluster.labels"]] <- singler.results$labels[match(seurat.obj[[]][["my.input.clusters"]], rownames(singler.results))]]


#custom heatmap
# Load necessary libraries
library(pheatmap)
library(RColorBrewer)

# Transpose your data and create a matrix
data_matrix <- matrix(
  c(0.6800488, 0.7615001, 0.7463142, 0.7249876, 0.6944787, 0.5765084, 0.7533016, 0.7734264, 0.7202637, 0.7170771, 0.7151528,
    0.6797350, 0.7504617, 0.7574219, 0.7242496, 0.6805648, 0.5647758, 0.7549216, 0.7778957, 0.7270300, 0.7209175, 0.7255079,
    0.6856680, 0.7439964, 0.7596284, 0.7260723, 0.6772450, 0.5726767, 0.7540155, 0.7718157, 0.7341739, 0.7120610, 0.7314577),
  nrow = 11, byrow = TRUE)

# Correcting the column names to match the number of rows in the transposed matrix
colnames(data_matrix) <- paste0("cluster_", 0:2)

# Row names to match the original columns
rownames(data_matrix) <- c("GW20_34_GE", "GW20_34_hypothalamus", "GW20_34_neocortex", "GW20_34_striatum", 
                           "GW20_34_thalamus", "GW20_claustrum", "GW20_GE", "GW20_hypothalamus", "GW20_neocortex", 
                           "GW20_striatum", "GW20_thalamus")

# Custom color palette
color_palette <- colorRampPalette(c("white", "red", "red1", "red2", "red3", "red4"))(100)
color_palette <- brewer.pal(n = 100, name = "Reds")
color_palette <- coolwarm(n=25)

# Find the column with the highest average score
highest_score_column <- which.max(colMeans(data_matrix))

new_pseudo1 <- t(new_pseudo1)
colnames(new_pseudo1) <- paste0("cluster ", 0:2)

# Generate annotations for rows and columns
annotation_col = data.frame(
  Structure = c("GW20_hypothalamus", "GW20_hypothalamus", "GW20_hypothalamus"))
rownames(annotation_col) = paste("cluster ", 0:2, sep = "")

# Specify colors
ann_colors = list(
  Structure = c("GW20_hypothalamus" = "orange")
)


# Creating the heatmap
pheatmap(new_pseudo1,
         color = color_palette,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         fontsize_row = 10,
         fontsize_col = 10,
         border_color = NA,
         legend_labels = c("Lower", "Higher"),
         legend_breaks = c(0,1)
) 


