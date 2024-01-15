
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Documents/Ruiqi_scRNAseq/10x_July2023/DS-Ngn2/filtered_feature_bc_matrix/")

# load expression matrix for editing
expression_matrix <- ReadMtx(
  mtx = "~/Documents/Ruiqi_scRNAseq/10x_July2023/DS-Ngn2/filtered_feature_bc_matrix/matrix.mtx.gz", features = "~/Documents/Ruiqi_scRNAseq/10x_July2023/DS-Ngn2/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "~/Documents/Ruiqi_scRNAseq/10x_July2023/DS-Ngn2/filtered_feature_bc_matrix/barcodes.tsv.gz"
)
expression_matrix <- ReadMtx(
  mtx = "~/Documents/Ruiqi_scRNAseq/wk4_filtered_feature_bc_matrix/matrix.mtx.gz", features = "~/Documents/Ruiqi_scRNAseq/wk4_filtered_feature_bc_matrix/features.tsv.gz",
  cells = "~/Documents/Ruiqi_scRNAseq/wk4_filtered_feature_bc_matrix/barcodes.tsv.gz"
)

#remove mm10 genes, select for GRCh38 only 
exp_mtx <-expression_matrix[stringr::str_detect(rownames(expression_matrix), "GRCh38"), ]
#change gene names to remove GRCh38
write.csv(rownames(exp_mtx), "~/Desktop/wk4_gene_names.csv")
genes <- read.csv("~/Desktop/wk4_gene_names.csv", header = FALSE)
rownames(exp_mtx) <- genes$V1

# Initialize the Seurat object with the raw (non-normalized data).
AD7 <- CreateSeuratObject(counts = exp_mtx, project = "AD7", min.cells = 3, min.features = 200)
AD7
AD3 <- CreateSeuratObject(counts = exp_mtx, project = "AD3", min.cells = 3, min.features = 200)
Ngn2 <- CreateSeuratObject(counts = exp_mtx, project = "Ngn2", min.cells = 3, min.features = 200)
Wk4 <- CreateSeuratObject(counts = exp_mtx, project = "wk4", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AD7[["percent.mt"]] <- PercentageFeatureSet(AD7, pattern = "^MT-")
AD3[["percent.mt"]] <- PercentageFeatureSet(AD3, pattern = "^MT-")
Ngn2[["percent.mt"]] <- PercentageFeatureSet(Ngn2, pattern = "^MT-")
Wk4[["percent.mt"]] <- PercentageFeatureSet(Wk4, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(AD3@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(wk4_vgat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(AD3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AD3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

AD7 <- subset(AD7, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
AD7 <- NormalizeData(AD7, normalization.method = "LogNormalize", scale.factor = 10000)

AD3 <- subset(AD3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
AD3_15 <- subset(AD3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)
Ngn2 <- subset(Ngn2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
test <- subset(wk4_vgat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < 10)

#merge samples
mature_INs <- merge(AD3, y = Wk4, add.cell.ids = c("AD3", "Wk4"), project = "mature_INs")
mature_INs
Wk4 <- NormalizeData(Wk4, normalization.method = "LogNormalize", scale.factor = 10000)

#variable features
Wk4 <- FindVariableFeatures(Wk4, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mature_INs), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Wk4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(Wk4)
Wk4 <- ScaleData(Wk4, features = all.genes)

#PCA
Wk4 <- RunPCA(Wk4, features = VariableFeatures(object = Wk4), dims = 1:20)
# Examine and visualize PCA results a few different ways
print(gcdata[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(gcdata, dims = 1:2, reduction = "pca")
DimPlot(gcdata, reduction = "pca")
DimHeatmap(gcdata, dims = 1, cells = 500, balanced = TRUE)

#determine dimensionality of dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
gcdata <- JackStraw(gcdata, num.replicate = 100)
gcdata <- ScoreJackStraw(gcdata, dims = 1:20)
JackStrawPlot(gcdata, dims = 1:20)

ElbowPlot(Wk4)

#cluster cells
Wk4 <- FindNeighbors(Wk4, dims = 1:20)
Wk4 <- FindClusters(Wk4, resolution = 0.2)
# Look at cluster IDs of the first 5 cells
head(Idents(gcdata), 5)

#run umap
Wk4 <- RunUMAP(Wk4, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Wk4, reduction = "umap", group.by = "seurat_clusters")



mature_INs_mito_stress_filt <- subset(x = mature_INs_mito_stress_filt, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)
Wk4 <- subset(x = Wk4, subset = SLC32A1 > 0)

#save object
saveRDS(AD3, file = "~/Desktop/AD3_new10x_mito20_20230926.rds")
save(mature_INs, file = "~/Desktop/mature_INs_mito20and10_20231003.RData")

#cell cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


AD3 <- CellCycleScoring(AD3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(AD3[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "orig.ident")

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
mature_INs_mito_stress_filt_slc32a1_or_gad_only <- RunPCA(mature_INs_mito_stress_filt_slc32a1_or_gad_only, features = c(s.genes, g2m.genes))
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only, reduction = "pca")

#Regress out cell cycle scores during data scaling
mature_INs_mito_stress_filt_slc32a1_or_gad_only <- ScaleData(mature_INs_mito_stress_filt_slc32a1_or_gad_only, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mature_INs_mito_stress_filt_slc32a1_or_gad_only))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed <- RunPCA(mature_INs_mito_stress_filt_slc32a1_or_gad_only, features = VariableFeatures(mature_INs_mito_stress_filt_slc32a1_or_gad_only), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed <- RunPCA(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed, features = c(s.genes, g2m.genes))
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed)
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed, reduction = "pca", group.by = "Phase")

save(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed, file = "~/Desktop/mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed_20231003.RData")

#set up for CSS integration
AD3 <- FindVariableFeatures(AD3, selection.method = "vst", nfeatures = 1000)
AD7 <- FindVariableFeatures(AD7, selection.method = "vst", nfeatures = 1000)
Ngn2 <- FindVariableFeatures(Ngn2, selection.method = "vst", nfeatures = 1000)
gcdata_mito_filtered <- FindVariableFeatures(gcdata_mito_filtered, selection.method = "vst", nfeatures = 1000)

write.csv(AD3@assays$RNA@var.features, file = "~/Desktop/AD3_varfeatures_1000.csv")
write.csv(AD7@assays$RNA@var.features, file = "~/Desktop/AD7_varfeatures_1000.csv")
write.csv(Ngn2@assays$RNA@var.features, file = "~/Desktop/Ngn2_varfeatures_1000.csv")
write.csv(gcdata_mito_filtered@assays$RNA@var.features, file = "~/Desktop/alltimepoints_new10x_varfeatures_1000.csv")

gcdata_mito_filtered <- FindClusters(gcdata_mito_filtered, resolution = 0.2)
DimPlot(gcdata_mito_filtered, reduction = "umap", group.by = "orig.ident")

#remove Ngn2
day2 <- subset(x = gcdata_CC_regressed_CSS_integrated, subset = timepoint == "48h_pD")


gcdata_mito_filtered_ADonly <- subset(x = gcdata_mito_filtered, subset = (orig.ident == "AD3" | orig.ident == "AD7" ))

OPRD1, OPRM1, OPRK1, OPR

FeaturePlot(object = wk4, reduction = "umap",features = c("OPRD1", "OPRM1", "OPRK1", "OPR"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)


ElbowPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed)
AD3_mito20_stress_filt_CC_and_ribo_regressed <- FindVariableFeatures(AD3_mito20_stress_filt_CC_and_ribo_regressed, nfeatures = 2000)
AD3_mito20_stress_filt_CC_and_ribo_regressed <- RunPCA(AD3_mito20_stress_filt_CC_and_ribo_regressed, features = VariableFeatures(object = AD3_mito20_stress_filt_CC_and_ribo_regressed), dims = 1:20)
#cluster cells
AD3_mito20_stress_filt_CC_and_ribo_regressed <- FindNeighbors(AD3_mito20_stress_filt_CC_and_ribo_regressed, dims = 1:20)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindClusters(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, resolution = .3)
# Look at cluster IDs of the first 5 cells
head(Idents(gcdata), 5)

#run umap
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- RunUMAP(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, reduction = "umap", group.by = "seurat_clusters", pt.size = .1)
save(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, file = "~/Desktop/gcdata_mito_10_wk2_15_filtered_stress_0.56_CC_regressed_20230817.RData")

g <- FindAllMarkers(AD3_mito20_stress_filt_CC_and_ribo_regressed, assay = "RNA")
g_up <- g[which(g$avg_log2FC > 0),]
g_up_sig <- g_up[which(g_up$p_val_adj < .05),]
View(g_up_sig)

write.csv(g, file = "~/Documents/Ruiqi_scRNAseq/new_10x_july2023/AD3_findallmarkers_DE_genes_6_clusters_20231026.csv")
write.csv(g_up, file = "~/Documents/Ruiqi_scRNAseq/new_10x_july2023/AD3_findallmarkers_DE_genes_6_clusters_up_only_20231026.csv")
write.csv(g_up_sig, file = "~/Documents/Ruiqi_scRNAseq/new_10x_july2023/AD3_findallmarkers_DE_genes_6_clusters_sig_up_only_20231026.csv")

FeaturePlot(object = AD3_filtered, features = c("SIX3", "SYT4", "EMX2", "SLC17A6", "SLC17A7", "LHX2", "LHX9", "GABRB2", "ONECUT2", "GRIA4", "ASCL1", "LHX1", "SST", "SOX2", "MAP2", "LMO3", "POU2F2", "FGF12", "LMO4", "TAC1", "RBFOX1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

test2 <- subset(gcdata_CC_regressed_CSS_integrated, subset = nFeature_RNA > 600)
DimPlot(test2, reduction = "umap_css", group.by = "timepoint")
FeatureScatter(gcdata_CC_regressed_CSS_integrated_stress_0.55_default, "nFeature_RNA", "percent.ribo", group.by = "timepoint", pt.size = 0.15)

gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- RunPCA(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, npcs = 30, verbose = FALSE)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- RunUMAP(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, reduction = "pca", dims = 1:30)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindNeighbors(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, reduction = "pca", dims = 1:30)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindClusters(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, resolution = 0.3)
DimPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, reduction = "pca", group.by = "timepoint")

FeaturePlot(object = Ngn2_filtered, features = c("GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN1", "GRIN3A", "GRIN3B", "MAP2", "SYN1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

wk4 <- subset(x = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "4wk_glia")

wk4_slc32a1_or_gad <- subset(x = wk4, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)
AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad <- subset(x = AD3_mito20_stress_filt_CC_and_ribo_regressed, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)

length(WhichCells(object = wk4, expression = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0))

save(AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad, file = "~/Desktop/AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad_10022023.RData")
save(wk4_slc32a1_or_gad, file = "~/Desktop/wk4_mito10_stress_filt_slc32a1_or_gad_10022023.RData")

#merge samples
mature_INs <- merge(AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad, y = wk4_slc32a1_or_gad, add.cell.ids = c("AD3", "wk4"), project = "mature_IN_integration")
#normalize
mature_INs <- NormalizeData(mature_INs, normalization.method = "LogNormalize", scale.factor = 10000)
#variable features
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindVariableFeatures(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, selection.method = "vst", nfeatures = 2000)
#scale data
all.genes <- rownames(mature_INs)
mature_INs <- ScaleData(mature_INs, features = all.genes)
#PCA
mature_INs <- RunPCA(mature_INs, features = VariableFeatures(object = mature_INs), dims = 1:20)
ElbowPlot(mature_INs)
#cluster cells
mature_INs <- FindNeighbors(mature_INs, dims = 1:20)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindClusters(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, resolution = 0.315)
#run umap
mature_INs <- RunUMAP(mature_INs, dims = 1:20)
DimPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, reduction = "umap", group.by = "seurat_clusters")
#save object
save(mature_INs, file = "~/Desktop/AD3_and_wk4_slc32a1_or_gad_only_10022023.RData")
saveRDS(mature_INs, file = "~/Desktop/AD3_and_wk4_slc32a1_or_gad_only_10022023.rds")


mature_INs$VIP_exprs <- NULL
object[['column.to.remove']] <- NULL


#set up for CSS integration
AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad <- FindVariableFeatures(AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad, selection.method = "vst", nfeatures = 1000)
wk4_slc32a1_or_gad <- FindVariableFeatures(wk4_slc32a1_or_gad, selection.method = "vst", nfeatures = 1000)
mature_INs <- FindVariableFeatures(mature_INs, selection.method = "vst", nfeatures = 1000)

write.csv(AD3_mito20_stress_filt_CC_regressed_slc32a1_or_gad@assays$RNA@var.features, file = "~/Desktop/AD3_slc32a1_or_gad_varfeatures_1000.csv")
write.csv(wk4_slc32a1_or_gad@assays$RNA@var.features, file = "~/Desktop/wk4_slc32a1_or_gad_varfeatures_1000.csv")
write.csv(mature_INs@assays$RNA@var.features, file = "~/Desktop/mature_INs_varfeatures_1000.csv")


FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("RPL15", "RPS13", "RPS11", "RPL10", "RPL27A", "RPS4X", "RPL34", "RPL34", "RPL32", "RPL41"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE) + DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters")

FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("ONECUT2", "MEF2C", "ARX", "TTC3", "RBFOX1", "SPOCK3", "MEIS2", "NCAM2", "SNHG6", "COX7C"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters")

AD3_mito20_stress_filt_CC_regressed <- PercentageFeatureSet(AD3_mito20_stress_filt_CC_regressed, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")

VlnPlot(AD3_mito20_stress_filt_CC_regressed, features = c("percent.mt", "percent.ribo"), ncol = 3)


#remove ribo genes
AD3_mito20_stress_filt_CC_regressed_ribo_removed <- AD3_mito20_stress_filt_CC_regressed[grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(AD3_mito20_stress_filt_CC_regressed), value = T, invert = T), ]

#regress out ribo genes
Ribo.small.subunit = toupper(c("RPSA", "rps2", "rps3", "rps3a", "rps4x", "rps4y", "rps5", "rps6", "rps7", "rps8", "rps9", "rps10", "rps11", "rps12", "rps13", "rps14", "rps15", "rps15a", "rps16", 'rps17', 'rps18', 'rps19', 'rps20', 'rps21', 'rps23', 'rps24',
                        'rps25', 'rps26', 'rps27', "rps27a", 'rps28', 'rps29', 'rps30'))

Ribo.large.subunit = toupper(c('rpl3','rpl4', 'rpl5', 'rpl6', 'rpl7', 'rpl7a', 'rpl8', 'rpl9','rpl10', 'rpl10a', 'rpl11','rpl12','rpl13', 'rpl13a','rpl14','rpl15','rpl16','rpl17','rpl18', 'rpl18a','rpl19','rpl21','rpl22','rpl23', 'rpl23a','rpl24','rpl26', 'rpl27',  'rpl27a', 'rpl28',
                       'rpl29', 'rpl30',  'rpl31','rpl32','rpl34','rpl35','rpl35a','rpl36', 'rpl36a', 'rpl37', 'rpl37a', 'rpl38', 'rpl39', 'rpl40', 'rpl41', 'rplp0', 'rplp1','rplp2'))



Ribo.genes = list(Ribo.small.subunit, Ribo.large.subunit)


mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed <- AddModuleScore(
  object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed,
  features = Ribo.genes,
  name = 'RiboScore',
  assay = "RNA"
)

mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed <- ScaleData(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed, vars.to.regress = c("RiboScore1", "RiboScore2"), features = rownames(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed))

AD3_mito20_stress_filt_CC_and_ribo_regressed <- RunPCA(AD3_mito20_stress_filt_CC_regressed, features = VariableFeatures(AD3_mito20_stress_filt_CC_regressed), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
AD3_mito20_stress_filt_CC_and_ribo_regressed <- RunPCA(AD3_mito20_stress_filt_CC_and_ribo_regressed, features = c(Ribo.small.subunit, Ribo.large.subunit))
DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed)
DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, reduction = "pca")

save(AD3_mito20_stress_filt_CC_and_ribo_regressed, file = "~/Desktop/AD3_mito20_stress_filtered_CC_and_ribo_regressed_20230928.RData")

#AD3 cluster specific DE feature plots
#cluster 0
FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("ONECUT2", "MEF2C", "ARX", "SDC2", "EGR1", "DLX6-AS1", "GABRA2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters", pt.size = .4)
#cluster 1
FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("COX7C", "UBL5", "COX17", "OST4", "NDUFA4", "ATP5ME", "SNHG9"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters", pt.size = .4)
#cluster 2
FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("LMO4", "MEIS2", "SPOCK3", "TAC1", "PTN", "RBFOX1", "NCAM2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters", pt.size = .4)
#cluster 3
FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("ZFHX3", "POU2F2", "ZFHX4", "SCGN", "CRNDE", "CRABP1", "DACH1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters", pt.size = .4)
#cluster 4
FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("HMX3", "BNC2", "RGMA", "ID4", "PLD5", "MAF", "ID2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters", pt.size = .4)
#cluster 5
FeaturePlot(object = AD3_mito20_stress_filt_CC_regressed, features = c("EMX2", "LHX5-AS1", "FSTL5", "PAX6", "TBR1", "SLC17A6", "FABP7"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)+ DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "seurat_clusters", pt.size = .4)



# Check for missing values (zero values)
missing_values <- rowMeans(mature_INs$RNA) == 0
table(missing_values)

# Filter the data to isolate genes with missing values (TRUE)
genes_with_missing_values <- mature_INs$RNA[missing_values, ]

# Remove missing value genes from data
mature_INs_genes_removed <- mature_INs[-genes_with_missing_values, ]



#create anndata file 
SaveH5Seurat(mature_INs, filename = "~/Desktop/mature_INs.h5Seurat")
Convert("~/Desktop/mature_INs.h5Seurat", dest = "h5ad")



library(SeuratData)
library(SeuratDisk)

###### WORKAROUND ######
# assigning the previous version of the `[[` function for the Assay class to the SeuratDisk package environment
"[[.Assay" <- function(x, i, ..., drop = FALSE) {
  if (missing(x = i)) {
    i <- colnames(x = slot(object = x, name = 'meta.features'))
  }
  data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
  if (drop) {
    data.return <- unlist(x = data.return, use.names = FALSE)
    names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
  }
  return(data.return)
}
environment(`[[.Assay`) <- asNamespace("SeuratObject")
rlang::env_unlock(asNamespace("SeuratDisk"))
assign("[[.Assay", `[[.Assay`, asNamespace("SeuratDisk"))
lockEnvironment(asNamespace("SeuratDisk"), bindings = TRUE)
rm(`[[.Assay`)
###### WORKAROUND ######

# example data
data("pbmc3k.final"); pbmc3k.final@images <- list()

# work
SaveH5Seurat(pbmc3k.final, "pbmc3k.final.h5Seurat", overwrite = T)
SaveH5Seurat(mature_INs, "~/Desktop/mature_INs.h5Seurat", overwrite = T)
#> Creating h5Seurat file for version 3.1.5.9900
#> Adding counts for RNA
#> Adding data for RNA
#> Adding scale.data for RNA
#> Adding variable features for RNA
#> Adding feature-level metadata for RNA
#> Adding cell embeddings for pca
#> Adding loadings for pca
#> No projected loadings for pca
#> Adding standard deviations for pca
#> Adding JackStraw information for pca
#> Adding cell embeddings for umap
#> No loadings for umap
#> No projected loadings for umap
#> No standard deviations for umap
#> No JackStraw data for umap

# This modification does not affect the Seurat and SeuratObject packages, so we can use the Seurat-v5 `[[` function
mature_INs[["RNA"]][[]]
#> [1] "counts"     "data"       "scale.data" 


test <- subset(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_regressed, subset = orig.ident == "wk4")
wk4_interneuron <- FindVariableFeatures(wk4_interneuron, nfeatures = 2000)
wk4_interneuron <- RunPCA(wk4_interneuron, features = VariableFeatures(object = wk4_interneuron), dims = 1:20)
#cluster cells
wk4_interneuron <- FindNeighbors(wk4_interneuron, dims = 1:20)
wk4_interneuron <- FindClusters(wk4_interneuron, resolution = .4)
wk4_interneuron <- RunUMAP(wk4_interneuron, dims = 1:20)
DimPlot(wk4_interneuron, reduction = "umap", group.by = "seurat_clusters", pt.size = .3)


wk4 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "4wk_glia")
wk4_IN_genes <- subset(wk4, subset = GAD2 > 0 | GAD1 > 0 | SLC32A1 > 0)
wk4_slc32a1_only2 <- subset(wk4_IN_genes, subset = SLC32A1 >0)

sum(GetAssayData(object = wk4_IN_genes, slot = "counts")["MAP2",] >0)
sum(GetAssayData(object = wk4_slc32a1_only, slot = "counts")["DCX",] >0)

g <- FindAllMarkers(wk4_map2_slc32a1_or_gad, assay = "RNA")
g_up <- g[which(g$avg_log2FC > 0),]
g_up_sig <- g_up[which(g_up$p_val_adj < .05),]
View(g_up_sig)

write.csv(g, file = "~/Documents/Ruiqi_scRNAseq/Week4_slc32a1_only_data_oct2023/wk4_map2_slc32a1_gad_only_findallmarkers_DE_genes_9_clusters_20231031.csv")
write.csv(g_up, file = "~/Documents/Ruiqi_scRNAseq/Week4_slc32a1_only_data_oct2023/wk4_map2_slc32a1_gad_only_findallmarkers_DE_genes_9_clusters_up_only_20231031.csv")
write.csv(g_up_sig, file = "~/Documents/Ruiqi_scRNAseq/Week4_slc32a1_only_data_oct2023/wk4_map2_slc32a1_gad_only_findallmarkers_DE_genes_9_clusters_sig_up_only_20231031.csv")
