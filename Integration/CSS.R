
devtools::install_github("quadbiolab/simspec")
library(simspec)
library(dplyr)

NormalizeData(gcdata_mito_filtered) %>%
  FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.9) %>%
  RunUMAP(dims = 1:20) %>%
  UMAPPlot(group.by = "seurat_clusters") 


NormalizeData(gcdata_mito_filtered) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 1) %>%
  RunUMAP(dims = 1:20) %>%
  UMAPPlot(group.by = "seurat_clusters") 

gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindNeighbors(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, dims = 1:20)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindClusters(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, resolution = 0.3) 
test <- RunUMAP(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, dims = 1:20)
test<- RunPCA(test, npcs = 20)
DimPlot(test, reduction = "umap", group.by = "timepoint")
UMAPPlot(gcdata_mito_filtered, group.by = "timepoint") 
gcdata_mito_filtered <- FindVariableFeatures(gcdata_mito_filtered, nfeatures = 2000)

load("~/Desktop/gcdata_mito_10and15_filtered_CC_regressed_20230518.RData")


genes <- read.csv("~/Desktop/union_top1000_var_genes_for_CSS_all_timepoints_20230820.csv", stringsAsFactors = FALSE)
genes <- read.csv("~/Desktop/union top 1000 genes 2023-05-17.csv", stringsAsFactors = FALSE)
genes <- read.csv("~/Desktop/union top 1000 genes with dualSMAD 2023-06-01.csv", stringsAsFactors = FALSE)
genes <- read.table("~/Desktop/union top 1000 genes with dualSMAD 2023-06-01.txt", stringsAsFactors = FALSE, header = TRUE)
genes <- read.csv("~/Desktop/AD3_and_wk4_union_1000_var_genes_for_CSS.csv", stringsAsFactors = FALSE)

gene <- genes$genes
gene <- genes$MARCH1
gene <- gcdata_mito_filtered@assays$RNA@var.features

#The calculation of CSS is implemented in the cluster_sim_spectrum function. 
#When applying to a Seurat v3 object, two fields are mandatory: object and label_tag. 
#object should be the provided Seurat object, the label_tag is the name of a column in the meta.data table that marks the sample/batch information for integration. 
#Here we calculate CSS to integrate organoids. It will take quite a while (approx 15-20 min).

mature_INs_CSS_integrated <- cluster_sim_spectrum(mature_INs, label_tag = "orig.ident", cluster_resolution = 0.4, var_genes = gene)
mature_INs_CSS_integrated <- RunUMAP(mature_INs_CSS_integrated, reduction="css", dims = 1:ncol(Embeddings(mature_INs_CSS_integrated,"css")), reduction.name="umap_css", reduction.key="UMACSS_")

save(mature_INs_CSS_integrated, file = "~/Desktop/AD3_wk4_slc32a1_or_gad_stress_filt_CC_regressed_CSS_integrated_1000_union_genes_res_0.4_20231002.RData")

mature_INs_CSS_integrated <- RunUMAP(mature_INs_CSS_integrated, reduction = "css", dims = 1:ncol(Embeddings(mature_INs_CSS_integrated, "css")))
mature_INs_CSS_integrated <- FindNeighbors(mature_INs_CSS_integrated, reduction = "css", dims = 1:ncol(Embeddings(mature_INs_CSS_integrated, "css")))
mature_INs_CSS_integrated <- FindClusters(mature_INs_CSS_integrated, resolution = 0.4)
UMAPPlot(mature_INs, reduction = "pca", group.by = "timepoint") + UMAPPlot(mature_INs, group.by = "timepoint") 
DimPlot(mature_INs, reduction = "umap_css", group.by = "orig.ident", cols = my_cols) + DimPlot(mature_INs, reduction = "umap_css", group.by = "seurat_clusters")
DimPlot(mature_INs_CSS_integrated, reduction = "umap_css", group.by = "seurat_clusters")

 my_cols <- c('48h_pD'='#E41A1C','7d_postDox'='#FF7F00','1wk_glia'='#E6AB02','2wk_glia'='#4DAF4A','3wk_glia'='#377EB8',
             '4wk_glia'='#984EA3','dualSMAD'='#A65628')

my_cols <- c('48h_pD'='#E41A1C','7d_postDox'='#FF7F00','1wk_glia'='#E6AB02','2wk_glia'='#4DAF4A','3wk_glia'='#377EB8',
             '4wk_glia'='#984EA3')

g1_untreat <- WhichCells(gcdata_CC_regressed_CSS_integrated, idents = "dualSMAD")
DimPlot(gcdata_CC_regressed_CSS_integrated, label=T, group.by="timepoint", cells.highlight= g1_untreat, cols.highlight = c("darkred", cols= "grey"))

DimPlot(gcdata_CC_regressed_CSS_integrated, reduction = "umap_css", group.by = "orig.ident", cols = my_cols) + DimPlot(gcdata_CC_regressed_CSS_integrated, reduction = "umap_css", group.by = "seurat_clusters")

        