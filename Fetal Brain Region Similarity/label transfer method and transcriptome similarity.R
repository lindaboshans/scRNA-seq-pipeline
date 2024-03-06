
PGP1_AD7_GW20_fetal_brain_seurat_integrated <- readRDS("~/Desktop/PGP1_AD7_GW20_fetal_brain_seurat_integrated_immune_anchors_20240229.rds")
AD7_map2_slc32a1 <- readRDS("~/Desktop/PGP1_WT2_AD7_map2_and_slc32a1_only_mito5_CC_regressed_20250228.rds")

PGP1_AD7_GW20_fetal_brain_seurat_integrated_no_AD7 <- subset(PGP1_AD7_GW20_fetal_brain_seurat_integrated, subset = dataset == "PGP1_AD7", invert = TRUE)


#Transcriptome similarity on cell cluster level using reference
avg_expr_ref <- sapply(sort(unique(PGP1_AD7_GW20_fetal_brain_seurat_integrated_no_AD7$structure_collapsed)), function(ct) rowMeans(PGP1_AD7_GW20_fetal_brain_seurat_integrated_no_AD7@assays$RNA$data[,which(PGP1_AD7_GW20_fetal_brain_seurat_integrated_no_AD7$structure_collapsed == ct)] ))
avg_expr_ds1 <- sapply(levels(AD7_map2_slc32a1$RNA_snn_res.0.4), function(ct) rowMeans(AD7_map2_slc32a1@assays$RNA$data[,which(AD7_map2_slc32a1$RNA_snn_res.0.4 == ct)]))

genes2cor <- intersect(VariableFeatures(PGP1_AD7_GW20_fetal_brain_seurat_integrated_no_AD7), rownames(AD7_map2_slc32a1))
corr2ref_cl <- cor(avg_expr_ds1[genes2cor,], avg_expr_ref[genes2cor,], method="spearman")

library(gplots)
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))


#transcitome similarity on single cell level
ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)
devtools::install_github("immunogenomics/presto")
library(presto)
ranked_expr_ds1 <- rank_matrix(AD7_map2_slc32a1@assays$RNA$data[genes2cor,])$X_ranked
install.packages("Matrix")
library(Matrix)

rank_matrix <- function (mat) 
{
  if (is.matrix(mat) | is.data.frame(mat)) {
    ranked_mat <- apply(mat, 2, rank)
  }
  else {
    df_mat <- Matrix::summary(mat)
    dfs_mat <- split(df_mat, df_mat$j)
    df_mat_ranked <- do.call(rbind, lapply(dfs_mat, function(df) {
      num_zeros <- nrow(mat) - nrow(df)
      ranks_nonzero <- rank(df$x)
      df$x <- ranks_nonzero + num_zeros - (1 + num_zeros)/2
      return(df)
    }))
    ranked_mat <- sparseMatrix(i = df_mat_ranked$i, j = df_mat_ranked$j, 
                               x = df_mat_ranked$x, dims = dim(mat), dimnames = dimnames(mat))
  }
  return(ranked_mat)
}
ranked_expr_ds1 <- rank_matrix(AD7_map2_slc32a1@assays$RNA$data[genes2cor,])

devtools::install_github("cysouw/qlcMatrix")
library(qlcMatrix)
corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
AD7_map2_slc32a1$celltype_maxcor <- ct_maxcor


plot1 <- UMAPPlot(AD7_map2_slc32a1)
plot2 <- UMAPPlot(AD7_map2_slc32a1, group.by="celltype_maxcor")
plot1 | plot2

corr2ref_scaled <- scale(t(corr2ref_cell))
corr2ref_sum2cl <- t(sapply(levels(AD7_map2_slc32a1@active.ident), function(cl)
  rowMeans(corr2ref_scaled[,which(AD7_map2_slc32a1@active.ident == cl)]) ))

heatmap.2(corr2ref_sum2cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))


#label transfer method
anchors <- FindTransferAnchors(reference = PGP1_AD7_GW20_fetal_brain_seurat_integrated_no_AD7, query = AD7_map2_slc32a1, dims = 1:30, npcs = 30)
predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$celltype, dims = 1:30)
seurat_DS1$celltype_transfer <- predictions$predicted.id

plot1 <- UMAPPlot(seurat_DS1, label=T)
plot2 <- UMAPPlot(seurat_DS1, group.by="celltype_transfer", label=T)
plot1 | plot2






#label transfer method 2

gcdata <- RunPCA(gcdata)
gcdata <- RunUMAP(gcdata, assay = "RNA", dims = 1:20)
gcdata <- FindNeighbors(gcdata, dims = 1:15)
gcdata <- FindClusters(gcdata, dims = 1:15, resolution = 0.5)
gcdata<- RunUMAP(gcdata, assay = "RNA", dims = 1:15)
DimPlot(gcdata, reduction = "umap", group.by = "seurat_clusters", pt.size = .25) 
g <- FindAllMarkers(gcdata, assay = "RNA")


avg_expr_ref <- sapply(sort(unique(PGP1_AD7_GW20_fetal_brain_seurat_integrated$structure_collapsed)), function(ct) rowMeans(PGP1_AD7_GW20_fetal_brain_seurat_integrated@assays$RNA$data[,which(PGP1_AD7_GW20_fetal_brain_seurat_integrated$structure_collapsed == ct)] ))
avg_expr_ref <- sapply(levels(stri_11pcw$Cell_type), function(ct) rowMeans(stri_11pcw@assays$RNA@data[,which(stri_11pcw$Cell_type == ct)] ))
avg_expr_ds1 <- sapply(levels(gcdata@active.ident), function(ct) rowMeans(gcdata@assays$RNA@data[,which(gcdata@active.ident == ct)]))
genes2cor <- intersect(VariableFeatures(stri_11pcw), rownames(gcdata))
genes2cor <- intersect(rownames(stri_11pcw), rownames(gcdata))
corr2ref_cl <- cor(avg_expr_ds1[genes2cor,], avg_expr_ref[genes2cor,], method="pearson")

library(gplots)
heatmap.2(corr2ref_cl, scale="none", trace="none", dendrogram = "none", key=T, keysize=1, key.par=list(mar=c(3,3.5,2,0)), , key.title = 'Correlation',
          key.xlab = 'Value', margins=c(14,16),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))
legend("right", title = "Tissues",legend=c("group1","group2"), 
       fill=c("red","green"), cex=0.8, box.lty=0)

#transcriptome similarity on cellular level
ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)
library(presto)
ranked_expr_ds1 <- rank_matrix(gcdata@assays$RNA@data[genes2cor,])$X_ranked
library(qlcMatrix)
corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
gcdata$celltype_maxcor <- ct_maxcor
plot1 <- UMAPPlot(gcdata, label=T)
plot2 <- UMAPPlot(gcdata, group.by="celltype_maxcor", label=F)
plot1 | plot2

corr2ref_scaled <- scale(t(corr2ref_cell))
corr2ref_sum2cl <- t(sapply(levels(gcdata@active.ident), function(cl)
  rowMeans(corr2ref_scaled[,which(gcdata@active.ident == cl)]) ))
heatmap.2(corr2ref_cl, scale="none", trace="none", key=F, keysize=0.5, margins=c(15,17),
          labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
          col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))

anchors <- FindTransferAnchors(reference = AdultNeuroExNac_V2, query = gcdata, dims = 1:30, npcs = 30)
predictions <- TransferData(anchorset = anchors, refdata = AdultNeuroExNac_V2$EmbryoAdultNuclei, dims = 1:30, k.weight = 2)
gcdata$celltype_transfer <- predictions$predicted.id

plot1 <- UMAPPlot(seurat_DS1, label=T)
plot2 <- UMAPPlot(gcdata, group.by="celltype_transfer", label=T)
plot1 | plot2

UMAPPlot(All.integrated, group.by="sampletype", label=T)



avg_expr_ref <- sapply(sort(unique(PGP1_AD7_GW20_fetal_brain_seurat_integrated$structure_collapsed)), function(ct) {
  # Extracting data for the current cell type
  data_matrix <- GetAssayData(PGP1_AD7_GW20_fetal_brain_seurat_integrated, assay = "RNA", slot = "data")[, WhichCells(PGP1_AD7_GW20_fetal_brain_seurat_integrated, expression = structure_collapsed == ct)]
  
  # Ensuring the result is a matrix with at least two dimensions
  if (is.vector(data_matrix)) {
    # If the result is a vector (only one cell in the category), convert it to a matrix
    data_matrix <- matrix(data_matrix, nrow = length(data_matrix), ncol = 1)
  }
  
  # Calculating row means, ensuring 'na.rm = TRUE' to handle any missing values
  rowMeans(data_matrix, na.rm = TRUE)
})

