BiocManager::install("clustifyr")
library(clustifyr)
week4_GW18_umap_corr_clustifyR_2000_ref_genes_RNA_assay

load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/GW20_brain_regions_integrated_slc32a1_only_Badhuri.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/week4_w_cell_anno_and_brain_corr.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Ruiqi_scRNAseq/RData_files/GW20_brain_regions_edited_slc32a1_only.RData")
GW2sai <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/seuset_sai.rds")
load("~/Desktop/gcdata_CC_regressed_CSS_integrated.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/GW18_GW20_fetal_regions_from_xuran_not_int.RData")

GW20_brain_regions_edited_slc32a1_only <- subset(GW20_brain_regions_integrated_slc32a1_only, subset = region.ident == "CGE" | region.ident == "MGE" | region.ident == "LGE" | region.ident == "PFCVZ", invert = TRUE)
save(GW20_brain_regions_edited_slc32a1_only, file = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/GW20_brain_regions_edited_slc32a1_only.RData")

GW18_20_fetal_regions <- subset(neuron.full.seu, subset = individual.y == "GW18_2" | individual.y == "GW18" | individual.y == "GW20" | individual.y == "GW20_31" | individual.y == "GW20_34" )
save(GW18_20_fetal_regions, file = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/GW18_GW20_fetal_regions.RData")
GW18_only <- subset(GW18_20_fetal_regions, subset = individual.y == "GW18_2" | individual.y == "GW18")
GW20_only <- subset(GW18_20_fetal_regions, subset = individual.y == "GW20" | individual.y == "GW20_31" | individual.y == "GW20_34")
fetal_brain_integrated <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW20" | dataset == "GW20_34")
br_regions_inhibitory <- subset(neuron.full.seu, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)


#make reference matrix
fetal_ages <- seurat_ref(
  seurat_object = neuron.full.seu,        # SeuratV3 object
  cluster_col = "individual.x")    # name of column in meta.data containing cell identities
br_regions_inhibitory <- FindVariableFeatures(br_regions_inhibitory, nfeatures = 2000)

strucutre, area, individual.y
g <- VariableFeatures(neuron.full.seu)

week4_slc32a1_only_vector <- clustify(
  input = wk4_map2_slc32a1_only, # matrix of normalized scRNA-seq counts (or SCE/Seurat object)
  cluster_col = "timepoint", # name of column in meta.data containing cell clusters
  ref_mat = fetal_ages, # matrix of RNA-seq expression data for each cell type
  seurat_out = FALSE,
  query_genes = g,
  compute_method = "spearman"
  #vec_out = TRUE #(to get vector to add to metadata column)
   # list of highly varible genes identified with Seurat
)

res2 <- cor_to_call(
  cor_mat = week4_slc32a1_only_vector,                  # matrix correlation coefficients
  cluster_col = "timepoint", # name of column in meta.data containing cell clusters
  threshold = 0.4
)

plot_cor_heatmap(cor_mat = week4_slc32a1_only_vector, col = c("blue", "white", "red")) 

plot5 <- FeaturePlot(week4_slc32a1_only,
                     features = c("GAD1", "GAD2"), label = TRUE, repel = TRUE) 

plot5

plot1 <- UMAPPlot(gcdata_seurat_integrated, group.by="seurat_clusters", label=T)
plot2 <- UMAPPlot(week4_slc32a1_only, group.by="type", label=F)
plot2 <- UMAPPlot(gc4, group.by="scType_w_CM_panDB_3", label=F)
plot2 <- UMAPPlot(GW20_brain_regions_edited_slc32a1_only, group.by="orig.ident", label=F)
plot1 | plot2

new_names <- read.csv("~/Desktop/test.csv", stringsAsFactors = FALSE, row.names = 1)
new <- read.delim("~/Desktop/test.txt", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
View(new)

rownames(neuron.full.seu@assays$RNA@meta.features) <- new$V1
rownames(neuron.full.seu@assays$RNA@counts) <- new$V1
rownames(neuron.full.seu@assays$RNA@data) <- new$V1

write.csv(week4_slc32a1_only, file = "~/Desktop/day7_fetal_brain_corr_matrix_4clusters_intr_fetal_brain.csv")
