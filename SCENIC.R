# Set up R environment
library(Seurat)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
library(Signac)
install.packages("tidyverse")
library(tidyverse)
library(data.table)
devtools::install_github("aertslab/SCopeLoomR")
library(SCopeLoomR)

load(file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_w_default_phate_slingshot_pseudotime.RData")
day2 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "48h_pD")
day7 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "7d_postDox")
wk1 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "1wk_glia")
wk2 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "2wk_glia")
wk3 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "3wk_glia")



# Extract raw count matrix [can also use normalized]
dgem <- wk4_map2_slc32a1_only@assays$RNA@counts
exp_mtx <- as.matrix(dgem)
dim(dgem)
head(colnames(dgem)) # columns => cell IDs
cell.info <- wk4_map2_slc32a1_only@meta.data

# Compute logical vector for every gene reflecting whether
# it has more than zero count per cell & and is detected in at 1 per thousand of cells in the global dataset
c <- floor(dim(dgem)[2]*.001); c  # c = 50 cells
nonzero <- dgem > 0L
keep_genes <- rowSums(as.matrix(nonzero)) >= c
filtered_dgem <- dgem[keep_genes, ]
dim(filtered_dgem)

# Further exclude mitochondrial genes
idx <- grep("^MT-", rownames(filtered_dgem))
rownames(filtered_dgem)[idx]
filtered_dgem <- filtered_dgem[-idx, ]
dim(filtered_dgem) # 20320 genes x 50931 cells



### OPTIONAL (necessary iif using pySCENIC for dimensionality analyses)

# Extract default embedding (e.g. UMAP or PCA coordinates)
default.umap <- Embeddings(wk4_map2_slc32a1_only, reduction = "umap")
default.umap.name <- "UMAP"


### Create loom file

file.name <- "~/Desktop/wk4_map2_slc32a1_only_mito10_stress_filt_CC_regressed_from_seurat_R_11062023.loom"
project.title <- "wk4_map2_slc32a1_only"
build_loom(
  file.name = file.name,
  dgem = filtered_dgem,
  title = project.title,
  genome = "hg38",
  default.embedding = default.umap,
  default.embedding.name = default.umap.name
)


### OPTIONAL

# Add further information to loom file (e.g. other embeddings)
loom <- open_loom(file.path = "~/Documents/Rachel_scRNAseq/SCENIC/SCENIC_visualizing_looms/gc_1wk_2wk_3wk_scenic_visualize.loom", mode = "r+")

loom <- open_loom(file.path = "~/Desktop/wk4_map2_slc32a1_only_mito10_stress_filt_CC_regressed_from_seurat_R_11062023.loom", mode = "r+")



other.umap <- Embeddings(day2, reduction = "PHATE")
add_embedding(loom = loom, embedding = other.umap, name = "PHATE")

pca.reduc <- Embeddings(wk4_map2_slc32a1_only, reduction = "pca")
add_embedding(loom = loom, embedding = pca.reduc, name = "PCA")

# Add Seurat cluster definition
add_seurat_clustering(loom = loom, seurat = wk4_map2_slc32a1_only, 
                      seurat.clustering.prefix = "RNA_snn_res.", 
                      default.clustering.resolution = "res.0.2")


### Add some custom clustering

#just making something up
kclust <- sctype(day2@meta.data$scType_w_CM_panDB_4)$cluster
add_clustering(loom = loom,
               group = "Custom clustering",
               name = "kmeans 3",
               clusters = kclust)

# Add annotation (categorical variable)
add_col_attr(
  loom=loom,
  key = "timepoint",
  value=as.character(wk3@meta.data$timepoint),
  as.annotation=T
)

# Add annotation (categorical variable)
add_col_attr(
  loom=loom,
  key = "seurat_5clusters",
  value=as.character(wk4_map2_slc32a1_only@meta.data$RNA_snn_res.0.4),
  as.annotation=T
)

# Add annotation (categorical variable)
add_col_attr(
  loom=loom,
  key = "slingshot_pseudotime",
  value=gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@meta.data$phate_slingshot_pseudotime1,
  as.annotation=F,
  as.metric = T
)

# Add annotation (categorical variable)
add_col_attr(
  loom=loom,
  key = "pseudotime_bins_final",
  value=as.character(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@meta.data$bins),
  as.annotation=T,
)


close_loom(loom)





### Define study-specific parameters

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
library(SCENIC)

out_dir <- "~/Documents/Ruiqi_scRNAseq/week4_map2_slc32a1_only_SCENIC_11062023/"
out_files <- list.files(out_dir, full.names = TRUE)

seurat_path <- "/Volumes/Linda_Mac_backup_2/gcdata_day7_final_20220829.RData"
analysis <- "week4_map2_slc32a1_only_SCENIC_11062023"





### Load pySCENIC output loom file (regulons)

# from pySCENIC working directory:

path2loom <- out_files[grep("pyscenic_wk4_map2_slc32a1_only_stress_filtered_aucell_output_11062023.loom", out_files)]
# using *filtered* loom to extract regulons
loom <- open_loom(path2loom, mode = "r")

# Read pySCENIC information from loom file
# WARNING: column attribute names may be different depending on the SCopeLoomR/pySCENIC version used!

# => regulons
regulonsMat <- get_regulons(loom,  # as incidence matrix
                            column.attr.name = "Regulons")

regulons_motif <- SCENIC::regulonsToGeneLists(regulonsMat) 

dim(regulonsMat)
dim(regulons_motif)

write.csv(regulonsMat,
          file = paste0(out_dir, "/", analysis, "_pySCENIC_regulons_incidence_matrix.csv"))

regulons <- SCENIC::regulonsToGeneLists(regulonsMat)  # can convert to gene list for processing with R
head(regulons)

regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
write.csv(regulonAUC@assays@data$AUC,
          file = paste0(out_dir, "/", analysis, "_pySCENIC_regulonsAUC.csv"))
embeddings <- get_embeddings(loom)
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- SCENIC::regulonsToGeneLists(regulons_incidMat)
regulonAucThresholds <- get_regulon_thresholds(loom)
write.csv(regulonAucThresholds,
          file = paste0(out_dir, "/", analysis, "_pySCENIC_regulonsAUCThresholds.csv"))
cellClusters <- get_clusterings(loom)
cellAnno_timepoint <- get_cell_annotation(loom, annotations.columns = "timepoint")
cellAnno_seurat_2clusters <- get_cell_annotation(loom, annotations.columns = "seurat_2clusters")
cellAnno_seurat_3clusters <- get_cell_annotation(loom, annotations.columns = "seurat_3clusters")
cellAnno_pseudotime <- get_cell_annotation(loom, annotations.columns = "slingshot_pseudotime")
cellAnno_bins <- get_cell_annotation(loom, annotations.columns = "pseudotime_bins_final")

cellAnno_timepoint <- as.data.frame(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$timepoint)
colnames(cellAnno_timepoint) <- 'timepoint'

close_loom(loom)


### Add desired data to original Seurat object (RNA only)

load("/Volumes/Linda_Mac_backup_2/gcdata_day7_final_20220829.RData")

# Extract underlying matrix from regulonsAUC object
AUCmat <- AUCell::getAUC(regulonAUC)
AUCmat[1:10, 1:10]
rownames(AUCmat)

# Add as assay to Seurat object
wk4_map2_slc32a1_only[['pyscenicAUC']] <- CreateAssayObject(data = AUCmat)
saveRDS(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, file = paste0(out_dir, "/", analysis, "_pyscenic2seurat.rds"))

save(wk4_map2_slc32a1_only, file = "~/Desktop/wk4_map2_slc32a1_only_w_pyscenicAUC_20231106.RData")

#exploring results in R

BiocManager::install(c("AUCell", "RcisTarget"))
devtools::install_github("aertslab/SCENIC") 

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# For some of the plots:
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

packageVersion("SCENIC")

setwd("~/Desktop/")
vsnDir <- "."

scenicLoomPath <- file.path("~/Desktop/all_timepoints_mito_10and15_07032023_SCENIC_integrated.loom")
scenicLoomPath <- file.path("~/Desktop/wk2_scenic_visualize.loom")
motifEnrichmentFile <- file.path("~/Documents/Ruiqi_scRNAseq/week4_interneuron_NeuronAndGabascore_only_SCENIC_11082023/regulons_wk4_interneuron_stress_filtered_11072023.csv")
file.exists(scenicLoomPath)
file.exists(motifEnrichmentFile)
list.files() #  What is already in the current work dir?

library(SCopeLoomR)
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat+1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
cellClusters <- get_clusterings(loom)
cellAnno <- get_cell_annotation(loom, annotations.columns = "resolution_0.1")
close_loom(loom)

# In case you need help:
?get_regulons

### To check whether it was read properly:
length(regulons);  head(names(regulons))
regulonAUC
length(regulonAucThresholds)
plot(embeddings$`Seurat_UMAP`)

motifEnrichment <- data.table::fread(motifEnrichmentFile, header=T, skip=1)[-3,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

write.csv(motifEnrichment,
          file = paste0(out_dir, "/", analysis, "_pySCENIC_motifEnrichment.csv"))

regulonAUC

#Regulators for known cell types or clusters
# To start from clusters/cell types from Scanpy: 
head(cellClusters)

bins <- as.data.frame(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins)
rownames(bins) <- colnames(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed)
levels(x = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins) <- c("", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins <- as.factor(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins)

selectedResolution <- "resolution_0.4" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(bins), bins[,"gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins"]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

cellClusters$bins <- as.numeric(bins$`gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins`)

#to use resolution annotation instead of default clusters
cellsPerCluster <- split(rownames(cellAnno_seurat_3clusters), cellAnno_seurat_3clusters[, 1]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

#to use resolution annotation instead of default clusters
cellsPerCluster <- split(rownames(bins), bins[, 1]) 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

#reorder columns so it goes 1 to 20
bins$`gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins` <-as.numeric(bins$`gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins`)
bins=bins[order(bins$`gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins`),]


# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))


#OR
d <- getAUC(regulonAUC)
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(data.frame(d[,cells])))

#OR
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     + function(cells) rowMeans(data.frame(regulon_data[,cells])))



#Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

#remove 0 column in bins
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[, -c(1)]
regulonActivity_byCellType <- regulonActivity_byCellType[, -c(1)]

#reorder columns
col_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
regulonActivity_byCellType <- regulonActivity_byCellType[, col_order]
head(regulonActivity_byCellType)
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType_Scaled[, col_order]
head(regulonActivity_byCellType_Scaled)



#plot:

options(repr.plot.width=8, repr.plot.height=10) # To set the figure size in Jupyter
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity", cluster_columns=FALSE,
                                   row_names_gp=grid::gpar(fontsize=2))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later

#regulon activity by celltype
regulonActivity_byCellType <- sapply(split(rownames(cell.info), cell.info$class_v8),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

#to see exact values
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "Clusters", "RelativeActivity")
topRegulators$Clusters <- factor(as.character(topRegulators$Clusters))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)
viewTable(topRegulators, options = list(pageLength = 10))

write.csv(regulonActivity_byCellType,"~/Documents/Ruiqi_scRNAseq/week4_map2_slc32a1_only_SCENIC_11062023/regulonActivity_by_5_clustsers.csv", row.names = TRUE)
write.csv(regulonActivity_byCellType_Scaled,"~/Documents/Ruiqi_scRNAseq/week4_map2_slc32a1_only_SCENIC_11062023/regulonActivity_scaled_by_5_clustsers.csv", row.names = TRUE)
write.csv(topRegulators,"~/Documents/Ruiqi_scRNAseq/week4_map2_slc32a1_only_SCENIC_11062023/topRegulators_by_5_clusters_11062023.csv", row.names = TRUE)
write.csv(topRegulators,"~/Documents/Rachel_scRNAseq/SCENIC/topRegulators_byCluster_wk2_2clusters.csv", row.names = TRUE)


#celltype specific regulators
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=bins[colnames(regulonAUC), 'gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins'])

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellAnno_seurat_3clusters[colnames(regulonAUC), 1])

selectedResolution <- "resolution_0.1" # select resolution
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellAnno2[colnames(regulonAUC), selectedResolution])


write.csv(rss,"~/Documents/Ruiqi_scRNAseq/week4_map2_slc32a1_only_SCENIC_11062023/regulon_specificity_score_RSS_by_5_clusters.csv", row.names = TRUE)


## Showing regulons and cell types with any RSS > 0.01 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot, height = 1500, width=450)
#Compare the average regulon Activity heatmaps and the RSS score: Find a few TFs that have a similar profile in both (e.g. specific for a cell type) and a few with different. Do you understand why? e.g. in the averaged heatmap there are regulons that are active across many cell types, while the RSS highlights the ones that are exclussive for each cell type.


options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = "2") # cluster ID

# List of embeddings available:
cat(names(embeddings), sep="\n")

# Overview of these embeddings (see below for details)
regulonsToPlot <- "DLX1_(+)"
options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2, ceiling(length(names(embeddings))/2)))
for (selectedEmbedding in names(embeddings))
  AUCell::AUCell_plotTSNE(embeddings[[selectedEmbedding]], exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .5, 
                          sub=selectedEmbedding)
selectedEmbedding <- embeddings[["Seurat_UMAP"]] # change if desired...


tfsToPlot <- c("POU2F2", "ATF3", "MEIS1") 
fsToPlot <- c("DLX1", "DLX2", "DLX6") 
regulonsToPlot <- unlist(lapply(tfsToPlot, function(x) grep(paste0("^", x,"_"), rownames(regulonAUC), value=TRUE)))

options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2,3))
# Plot expression:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log[tfsToPlot,], plots=c("Expression"), cex = .4)
# Plot regulon activity:
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, regulonAUC[regulonsToPlot,], plots=c("AUC"), cex = .4)

#plot all regulons at once
# Subset some cells to make the plots faster:
nCells <- 3000
set.seed(123)
cellsSelected <- sample(colnames(regulonAUC), nCells) 
regulonAUC_subset <- regulonAUC[regulonOrder, which(colnames(regulonAUC) %in% cellsSelected)]
dim(regulonAUC_subset)
selectedEmbedding_subset <- selectedEmbedding[colnames(regulonAUC_subset), ]

# Save AUC as PDF:
pdf("wk2_RegulonActivity_ALL.pdf", width=20, height=15)
par(mfrow=c(4,6))

AUCell::AUCell_plotTSNE(selectedEmbedding_subset, cellsAUC=regulonAUC_subset, plots="AUC")
dev.off()

#binarized regulon activity
head(as.data.frame(regulonAucThresholds))

regulonAucThresholds <- setNames(names(regulonAucThresholds), regulonAucThresholds)

regulonsToPlot <- "DLX2_(+)"
options(repr.plot.width=10, repr.plot.height=4) # To set the figure size in Jupyter
par(mfrow=c(1,3))
AUCell::AUCell_plotTSNE(selectedEmbedding, exprMat_log, 
                        regulonAUC[regulonsToPlot,], thresholds = regulonAucThresholds[regulonsToPlot],
                        plots=c("AUC", "histogram", "binary"), cex = .5)

# This function will be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}

binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)
binaryRegulonActivity[1:5,1:3]

# Subset some cells to make the plots faster:
nCells <- 1000
set.seed(123)
cellsSelected <- sample(colnames(regulonAUC), nCells) 
binAct_subset <- binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% cellsSelected)]
dim(binAct_subset)

options(repr.plot.width=12, repr.plot.height=10) # To set the figure size in Jupyter
# binAct_subset <- binAct_subset[regulonOrder,]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Binarized activity", 
                        col = c("white","pink","red"),
                        cluster_rows = TRUE, cluster_columns=FALSE,
                        show_column_names = TRUE,
                        row_names_gp=grid::gpar(fontsize=6)) # row font size

minPerc <- .75
cellInfo_binarizedCells <- cellAnno_seurat_5clusters[which(rownames(cellAnno_seurat_5clusters)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
cellInfo_binarizedCells <- bins[which(rownames(bins)%in% colnames(binaryRegulonActivity)),, drop=FALSE]

regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seurat_5clusters),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))

#reorder columns
col_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
regulonActivity_byCellType_Binarized <- regulonActivity_byCellType_Binarized[, col_order]
head(regulonActivity_byCellType_Binarized)

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
dim(binaryActPerc_subset)
head(binaryActPerc_subset)

#remove 0 column 
binaryActPerc_subset <- binaryActPerc_subset[, -c(1)]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"), row_names_gp = grid::gpar(fontsize = 3))


ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Binarized activity", 
                        col = c("white", "pink", "red"),
                        cluster_rows = TRUE, cluster_columns=FALSE,
                        show_column_names = TRUE,
                        row_names_gp=grid::gpar(fontsize=6)) # row font size


AD3_mito20_stress_filt_CC_and_ribo_regressed[['AUCBinary']] <- CreateAssayObject(data = binaryRegulonActivity)


save(AD3_mito20_stress_filt_CC_and_ribo_regressed, file = "~/Desktop/AD3_mito20_stress_filt_CC_and_ribo_regressed_w_pyschenicAUC_and_AUCBinary_10022023.RData")


write.csv(regulonActivity_byCellType_Binarized,"~/Documents/Ruiqi_scRNAseq/alltimepoints_post_stress_SCENIC_09082023/by_pseuodtime_bins/newer_binary_heatmaps_10262023/regulonActivity_by_bins_Binarized_alltimepoints.csv", row.names = TRUE)
write.csv(binaryActPerc_subset,"~/Documents/Ruiqi_scRNAseq/alltimepoints_post_stress_SCENIC_09082023/regulonActivity_by_bins_Binarized_alltimepoints_minpercent_.75.csv", row.names = TRUE)


gene <- "GAD1"
names(regulons)[which(sapply(regulons, function(x) "GAD1" %in% x))]

dim(regulons_incidMat)

genes <- c("GAD1", "GAD2") 
incidMat_subset <- regulons_incidMat[,genes]
incidMat_subset <- incidMat_subset[rowSums(incidMat_subset)>0,]

write.csv(incidMat_subset, "~/Documents/Ruiqi_scRNAseq/alltimepoints_post_stress_SCENIC_09082023/GAD1_GAD2_containing_regulons.csv")

write.csv(binaryRegulonActivity,"~/Documents/Rachel_scRNAseq/SCENIC/regulonActivity_byCluster_48hr.csv", row.names = FALSE)
write.csv(regulonActivity_byCellType,"~/Documents/Rachel_scRNAseq/SCENIC/regulonActivity_byCluster_48hr.csv", row.names = FALSE)





df <- read.csv("~/Documents/Ruiqi_scRNAseq/week4_interneuron_NeuronAndGabascore_only_SCENIC_11082023/wk4_interneuron_regulonActivity_by_3_clustsers.csv", sep = ",", header = TRUE, row.names = 1)
regAct <- as.matrix(df)
mean <- mean(regAct)
sd <- sd(regAct)
> z_cutoff <- qnorm(0.90)

actual_cutoff <- mean + (z_cutoff * sd)


# Specify the z-score cutoff value
cutoff_value <- 0.1018992

# Filter the data frame based on the cutoff value for each column
filtered_df <- regAct
filtered_df[df < cutoff_value] <- NA

write.csv(filtered_df, file = "~/Documents/Ruiqi_scRNAseq/week4_interneuron_NeuronAndGabascore_only_SCENIC_11082023/wk4_interneuronregulonActivity_scaled_by_3_clustsers_top10_percent_regulons.csv")


#paint umap with regulon activity scores

# Assuming 'seurat_object' is your Seurat object and 'regulon_activity_scores' is a matrix or data frame
# with regulon activities (rows are cells, columns are regulons)
regulons <- read.csv("~/Documents/Ruiqi_scRNAseq/alltimepoints_post_stress_SCENIC_09082023/by_pseuodtime_bins/regulonActivity_scaled_by_pseuodotime_bins.csv", header = TRUE)
head(regulons)


# Add regulon activity scores to Seurat metadata
seurat_object[["regulon_activity"]] <- CreateAssayObject(counts = regulon_activity_scores)

# Choose a specific regulon to plot
regulon_to_plot <- "YourRegulonName"

# Normalize and scale the data if not already done
seurat_object <- NormalizeData(seurat_object, assay = "regulon_activity")
seurat_object <- ScaleData(seurat_object, assay = "regulon_activity")

# Run UMAP if not already done
seurat_object <- RunUMAP(seurat_object, dims = 1:10, assay = "regulon_activity")

# Plot UMAP with regulon activity
DimPlot(seurat_object, reduction = "umap", group.by = regulon_to_plot)

