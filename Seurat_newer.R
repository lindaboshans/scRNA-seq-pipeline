#Seurat
#instal libraries

install.packages('Seurat')

library(Seurat)
library(dplyr)
library(Matrix)
#load data
load(file = "/Volumes/Linda_Mac_backup_2/gcdata_48hr_final_20220830.RData")
load(file = "/Users/lindaseong/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/gcdata_day7_final_20220829.RData")
load(file = "~/Documents/Rachel_scRNAseq/gcdata_day7_final_20220829.RData")
load(file = "/Volumes/Linda_Mac_backup_2/Rachel_scRNAseq/gcdata_LB.RData")
read_rds("/Volumes/Linda_Mac_backup_2/Rachel_scRNAseq/seuset_sai.rds")
load(file = "/Users/lindalee/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_LB.RData")
load(file = "/Users/lindalee/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/dualSMADi_gcdata.rda")


# set identity 
Idents(gc_week4_neuron) <- 'GAD1.GAD2.SLC32A1'
Idents(gcdata_mito_filtered) <- 'timepoint'
Idents(gcdata_CC_regressed_CSS_integrated) <- 'timepoint'
gcdata_CC_regressed_CSS_integrated$timepoint <- factor(gcdata_CC_regressed_CSS_integrated$timepoint, levels = c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia', 'dualSMAD'))


#rename idents
gcdata <- RenameIdents(gcdata, 'MAP2-TUBB3-DCX-' = 'non-neuronal', 'MAP2-TUBB3-DCX+' = 'neuron', 'MAP2+TUBB3-DCX-' = 'neuron', 'MAP2+TUBB3-DCX+' = 'neuron', 'MAP2-TUBB3+DCX-' = 'neuron', 'MAP2+TUBB3+DCX-' = 'neuron', 'MAP2-TUBB3+DCX+' = 'neuron', 'MAP2+TUBB3+DCX+' = 'neuron')
gc_week4_neuron <- RenameIdents(gc_week4_neuron, 'GAD1+GAD2+SLC32A1+' = 'GABAergic', 'GAD1-GAD2+SLC32A1-' = 'GABAergic', 'GAD1-GAD2-SLC32A1+' = 'GABAergic', 'GAD1+GAD2-SLC32A1-' = 'GABAergic', 'GAD1+GAD2-SLC32A1+' = 'GABAergic', 'GAD1+GAD2+SLC32A1-' = 'GABAergic', 'GAD1-GAD2+SLC32A1+' = 'GABAergic', 'GAD1-GAD2-SLC32A1-' = 'non-GABA neuron')

#stash idents as new column
gc_week4_neuron[["GABA.ident"]] <- Idents(object = gc_week4_neuron)

#plot some initial genes
DotPlot(gcdata, features = c("ATF4", "ATF6", "XBP1", "HSPA5", "SLC12A5", "SLC12A2"), group.by = "idents.tp") + RotatedAxis() & scale_colour_gradientn(colours =rev(brewer.pal(n = 11, name = "RdBu")))

#scale data
all.genes <- rownames(gcdata)
gcdata <- ScaleData(gcdata, features = all.genes)

#make umap (all were done with dims 1:20)
gcdata <- RunUMAP(object = gcdata, assay = "RNA", dims = 1:20)                        
DimPlot(gcdata, dims = c(1,2), reduction = "umap", group.by = "seurat_clusters", pt.size = 1) 

    
 # To subset and remove single cluster and keep the remaining clusters for new analysis
sub_obj <- subset(x = gcdata$timepoint, idents = 1, invert = TRUE)
# If you have renamed the clusters be sure to provide their names in quotes in the function
gc <- subset(x = gcdata, idents = "interneuron.4wk_glia")
#subset for 2 conditions
gc_week4_neuron <- subset(x = gcdata, subset = (idents.tp == "neuron_4wk_glia" | idents.tp == "non-neuronal_4wk_glia"))
    
#subset to get proportion of SLC32A1+ cells out of "neurons" at 4 weeks.
gc <- subset(gc_week4_neuron, subset = celltype_lb == "neuron")
gc_week4_slc32a1_only <- subset(gc_week4_neuron, subset = GABA.ident == "GABAergic")
gc_week4_slc32a1_only <- subset(gc_week4_neuron, subset = (GAD1.GAD2.SLC32A1 == "GAD1-GAD2-SLC32A1+" | GAD1.GAD2.SLC32A1 == "GAD1+GAD2-SLC32A1+" | GAD1.GAD2.SLC32A1 == "GAD1-GAD2+SLC32A1+" | GAD1.GAD2.SLC32A1 == "GAD1+GAD2+SLC32A1+"))
    
#get % of cells expressing gene
sum(GetAssayData(object = gc_week4_neuron, slot = "counts")["SLC32A1",] >=1)
    
#change resolution of cluster
gcdata <- FindClusters(gcdata, resolution = 0.6)
    
#Dotplots of key genes 
TFs <- c("CHAT", " KISS1", " TH", " TRH", " AVP", " CCK", " VIP", " SIX6", " LHX1", " CALCA", " SOSTDC1", " TAC3", " TAC1", " NPW", " BDNF", " AMIGO2", " CHRM2", " SLC17A6", " PMAIP", " CRH", " LHX8", " POU3F3", " PNOC", " NR2F2", " IGSF1", " NRGN", " RELN", " NPY", " PTHLH", " PROK2", " PDYN", " DRD1", " PENK", " CXCL14", " ARPP21", " ETV1", " MOXD1", " LYPD6", " SST", " MYLK", " QRFP", " HCRT", " AGRP", " GNRH", " GRP", " GAL", " NMU", " NMB", " SCGN", " NTS", " GHRH", " ONECUT2", " MEIS2", " ISL1", " FOXP1", " FOXP2", " PBX1", " PBX3", " SP9", " SOX4", " SOX11")
DotPlot(wk4_slc32a1_only, features = c("CHAT", "KISS1", "TH", "TRH", "AVP", "CCK", "VIP", "SIX6", "LHX1", "CALCA", "SOSTDC1", "TAC3", "TAC1", "NPW", "BDNF", "AMIGO2", "CHRM2", "SLC17A6", "PMAIP", "CRH", "LHX8", "POU3F3", "PNOC", "NR2F2", "IGSF1", "NRGN", "RELN", "NPY", "PTHLH", "PROK2", "PDYN", "DRD1", "PENK", "CXCL14", "ARPP21", "ETV1", " MOXD1", " LYPD6", "SST", "MYLK", "QRFP", "HCRT", "AGRP", "GNRH", "GRP", "GAL", "NMU", "NMB", "SCGN", "NTS", "GHRH", "ONECUT2", "MEIS2", "ISL1", "FOXP1", "FOXP2", "PBX1", "PBX3", "SP9", "SOX4", "SOX11"))
    
DoHeatmap(gc4, features = c("CHAT", "KISS1", "TH", "TRH", "AVP", "CCK", "VIP", "SIX6", "LHX1", "CALCA", "SOSTDC1", "TAC3", "TAC1", "NPW", "BDNF", "AMIGO2", "CHRM2", "SLC17A6", "PMAIP", "CRH", "LHX8", "POU3F3", "PNOC", "NR2F2", "IGSF1", "NRGN", "RELN", "NPY", "PTHLH", "PROK2", "PDYN", "DRD1", "PENK", "CXCL14", "ARPP21", "ETV1", " MOXD1", " LYPD6", "SST", "MYLK", "QRFP", "HCRT", "AGRP", "GNRH", "GRP", "GAL", "NMU", "NMB", "SCGN", "NTS", "GHRH", "ONECUT2", "MEIS2", "ISL1", "FOXP1", "FOXP2", "PBX1", "PBX3", "SP9", "SOX4", "SOX11"), group.by = "seurat_clusters")

#resolution 1.2 
DotPlot(gc, features = c("TRH","VIP", "TAC3", "TAC1", "BDNF", "AMIGO2", "PMAIP", "IGSF1", "NRGN", "RELN", "NPY", "ARPP21", "ETV1", " MOXD1", " LYPD6", "SST", "GNRH", "GRP", "GAL", "NMB", "SCGN", "GHRH", "ONECUT2", "MEIS2", "ISL1"), group.by = "seurat_clusters") & theme(axis.text.x=element_text(angle=45, hjust=1, size=9)) & scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu")))
   
#hypothalamus genes 0.6
DotPlot(gc, features = c("GAL", "CALB1", "ONECUT2", "ARPP21", 	"NPY",	"SST", "CALY", "IGSF1", "AMIGO2",	"NRGN",	"VIP",	"GRP",	"SCGN",	"GHRH",	"RELN",	"TAC3",	"CALB2", "TAC1", "TRH", "MEIS2"), group.by = "seurat_clusters") & theme(axis.text.x=element_text(angle=45, hjust=1, size=9)) & scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu")))
    
DotPlot(wk4_slc32a1_only, features = c("GAL",	"SST", "TAC1",	"TRH",	"NPY",	"ARPP21",	"IGSF1",	"AMIGO2",	"NRGN",	"VIP",	"GRP",	"SCGN",	"GHRH",	"RELN",	"TAC3",	"CARBP1",	"CALY",	"SYT1",	"SYT4",	"CALB1",	"CALB2",	"SEC11C",	"SOX4",	"MEIS2",	"PCSK1N",	"GAP43"), group.by = "seurat_clusters") & theme(axis.text.x=element_text(angle=45, hjust=1, size=9)) & scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu")))
    
#Find cells expressing combo of genes and add as column to metadata
cells = colnames(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed)
WhichCells(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, cells = cells, expression = SLC32A1 > 1 & GAD1 < 1 & GAD2 < 1)
    
Neuron <- WhichCells(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, expression = MAP2 > 0)
    
gcdata$MAP2 <- ifelse(colnames(gcdata) %in% Neuron_lb, "MAP2+", "MAP2-")
    
#remove unwanted columns 
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed[['hash.ID ']] <- NULL

    
#module score for set of genes 
genes.for.scoring <- list(c("MAP2", "TUBB3", "DCX"))
genes.for.scoring <- list(c("GAD1", "GAD2", "SLC32A1"))
#calculate module score
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- AddModuleScore(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, features = genes.for.scoring, ctrl = 100, name = "neuron.score.lb", random.seed = 1)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- AddModuleScore(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, features = genes.for.scoring, ctrl = 100, name = "GABA.score.lb", random.seed =1)
    
#annotate module score column 
test <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@meta.data %>% mutate(neuronal.ident = if_else (neuron.score.lb1 >0, "neuron", "non-neuron"))
test <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@meta.data %>% mutate(GABA.ident = if_else (GABA.score.lb1 >0, "interneuron", "non-GABAergic"))
    
#add column to object
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed<- AddMetaData(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, metadata = test$neuronal.ident, col.name = 'neuron_lb')
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- AddMetaData(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, metadata = test$GABA.ident, col.name = 'GABA.ident_lb')
    
save(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, file = "~/Desktop/RData files/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_phate_slingshot_pseudotime_pyschenicAUC_AUCBinary_gaba_scores_11072023.RData")
save(wk4, file = "~/Desktop/RData files/wk4_allcells_with_gaba_and_neuron_module_scores_11072023.RData")
save(wk4_interneuron, file = "~/Desktop/RData files/wk4_interneurons_with_positive_gaba_and_neuron_module_scores_11072023.RData")
    
    #plot gene expression within clusters and change color
    FeaturePlot(object = gc, features = "SLC17A7", cols = c("gray", "red"))
    FeaturePlot(object = gc, features = "GAD1", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "GAD2", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "SLC17A7", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "SLC17A6", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
   FeaturePlot(object = gc, features = "TPH1", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "TPH2", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "TH", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "SLC6A5", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "CHAT", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "SLC18A3", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "SLC18A2", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    FeaturePlot(object = gc, features = "SLC32A1", slot = "counts", cols = c("gray", "red"), pt.size = 0.5, order = TRUE)
    
    
    
#phateR using seurat
reticulate::py_install("phate", pip=TRUE)
devtools::install_github("KrishnaswamyLab/phateR")
    
install.packages("phateR")
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
if (!require(Rmagic)) install.packages("Rmagic")
library(Rmagic)
    
devtools::install_github("scottgigante/seurat", ref="patch/add-PHATE-again")
    
### handles the Seurat to phate conversion
    
### first, grab the input required for phate (here we are using the normalized data stored in Seurat
seurat_data <- as.data.frame(dlx6pos@assays$RNA@data)
    
## reshape for input into PHATE
phate_data_input <- t(seurat_data)
    
## run phate 
phate_output <- phate(phate_data_input) 
    
## quick sanity check of the phate embeddings
ggplot(phate_output, aes(x=PHATE1, y=PHATE2)) +
    geom_point()
    
## stash the embeddings back into the seurat object as a dimension reduction object
seurat.object[["PHATE"]] <- CreateDimReducObject(embeddings = phate_output$embedding, key = "PHATE_", assay = DefaultAssay(seurat.object))
    
## plot using seurat's default tools
DimPlot(seurat.object , reduction = "PHATE") + ggtitle(label = "PHATE")
#plot with ggplot
plot <- ggplot(seuset@reductions$phate) + geom_point(aes(PHATE1, PHATE2, color=gcdata$timepoint), size= 0.5) + 
    labs(color="timepoint") + scale_color_viridis(discrete=TRUE) + theme_bw()
    

    #Data Integration using anchors
    gc48hr <- SplitObject(gcdata, split.by = "timepoint")
    gc48hr <- lapply(gc48hr, SCTransform)
    features  <- SelectIntegrationFeatures(gc48hr, nfeatures = 3000)
    gc48hr <- PrepSCTIntegration(gc48hr, anchor.features = features)
    anchors <- FindIntegrationAnchors(gc48hr, normalization.method = "SCT", anchor.features = features)
    gc48hr <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                                     k.weight = 46)
    gc48hr <- RunPCA(gc48hr)
    gc48hr <- RunUMAP(gc48hr, reduction = "pca", dims = 1:20)
    gc48hr <- FindNeighbors(gc48hr, dims = 1:20)
    gc48hr <- FindClusters(gc48hr, dims = 1:20, resolution = 0.2)
    gc48hr<- RunUMAP(gc48hr, assay = "RNA", dims = 1:20)
    DimPlot(gc4, reduction = "umap", group.by = "seurat_clusters", pt.size = .25)
 

#subcluster using harmony
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("harmony")

library(harmony)
gc_7d <- Seurat::NormalizeData(gc_7d, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = gc_7d@var.genes, npcs = 20, verbose = FALSE)

DimPlot(object = gc_7d, reduction = "pca", pt.size = .1, group.by = "seurat_clusters")


gc_7d <- RunHarmony(gc_7d, "timepoint", ndims = 1:25)
DimPlot(object = gc_7d, reduction = "harmony", pt.size = .1, group.by = "seurat_clusters")

gc_7d <- gc_7d %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()
DimPlot(gc_7d, reduction = "harmony", group.by = "seurat_clusters", pt.size = .1)

gc_7d <- RunHarmony(gcdata, "timepoint")
gc_7d <- RunUMAP(gc_7d, reduction = "harmony", dims = 1:20)





neuron.markers <- FindAllMarkers(gcdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
neuron.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

neuron.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(gcdata, features = top10$gene) + NoLegend()


#pull cells expression certain genes for DE
Idents(gcdata, WhichCells(object = gcdata, expression = ARX > 0, slot = 'data')) <- 'ARX.pos'
Idents(gcdata, WhichCells(object = gcdata, expression = GATA2 > 0, slot = 'data')) <- 'GATA2.pos'
Idents(gcdata, WhichCells(object = gcdata, expression = POU5F1 > 0, slot = 'data')) <- 'POU5F1.pos'
Idents(gcdata, WhichCells(object = gcdata, expression = `DLX6-AS1` > 0  & ARX <=0, slot = 'data')) <- 'DLX6as1.pos_ARX.neg'

#find differential expression between one cell identity and all other clusters
DE_genes_ARX_all <- FindMarkers(gcdata, ident.1 = 'ARX.pos', ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod', min.pct = 0.1, random.seed =1)
View(DE_genes_ARX_all)

# Combine markers with gene descriptions 
ann_comb_markers <- inner_join(x = DE_genes_ARX_all, 
                               y = annotations[, c("gene_name", "description")],
                               by = c("gene" = "gene_name")) %>%
  unique()
# Rearrange the columns to be more intuitive
ann_comb_markers <- ann_comb_markers[ , c(6, 7, 2:4, 1, 5,8)]

# Order the rows by p-adjusted values
ann_comb_markers <- ann_comb_markers %>%
  dplyr::arrange(cluster, p_val_adj)

View(ann_comb_markers)

# Save markers to file
write.csv(DE_genes_ARX_all, 
          file = "~/Documents/Rachel's scRNAseq/ARX_pos_cells_vs_all_others_DE.csv", 
          quote = FALSE, 
          row.names = TRUE)

#find differential expression between one cell identity and all other clusters
DE_genes_DLX6as1_all <- FindMarkers(gcdata, ident.1 = 'DLX6as1.pos_ARX.neg', ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod', min.pct = 0.1, random.seed =1)
View(DE_genes_DLX6as1_all)

# Order the rows by p-adjusted values
DE_genes_DLX6as1_all <- DE_genes_DLX6as1_all %>%
  dplyr::arrange(p_val_adj)

# Save markers to file
write.csv(DE_genes_DLX6as1_all, 
          file = "~/Documents/Rachel's scRNAseq/DLX6as1_pos_ARX_neg_cells_vs_all_others_DE.csv", 
          quote = FALSE, 
          row.names = TRUE)

#find differential expression between one cell identity and all other clusters
DE_genes_GATA2_all <- FindMarkers(gcdata, ident.1 = 'GATA2.pos', ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod', min.pct = 0.1, random.seed =1)
View(DE_genes_GATA2_all)

# Save markers to file
write.csv(DE_genes_GATA2_all, 
          file = "~/Documents/Rachel's scRNAseq/GATA2_pos_cells_vs_all_others_DE.csv", 
          quote = FALSE, 
          row.names = TRUE) 

#find differential expression between one cell identity and all other clusters
DE_genes_POU5F1_all2 <- FindMarkers(gcdata, ident.1 = 'POU5F1.pos', ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod', min.pct = 0.1, random.seed =1)
View(DE_genes_POU5F1_all2)

# Save markers to file
write.csv(DE_genes_POU5F1_all, 
          file = "~/Documents/Rachel's scRNAseq/POU5F1_pos_cells_vs_all_others_DE.csv", 
          quote = FALSE, 
          row.names = TRUE) 

#find differential expression between one cell identity and all other clusters
DE_genes_cluster_0vs1 <- FindMarkers(gcdata, ident.1 = 0, ident.2 = 1, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod', min.pct = 0.1, random.seed =1)
View(DE_genes_cluster_0vs1)

# Save markers to file
write.csv(DE_genes_cluster_0vs1, 
          file = "~/Documents/Rachel's scRNAseq/48hr_DE_genes_cluster_0vs1_res.08.csv", 
          quote = FALSE, 
          row.names = TRUE)


#find differential expression between one cell identity and all other clusters
DE_genes_cluster_4vsall <- FindMarkers(objct = gcdata, ident.1 = 4, ident.2 = NULL, only.pos = TRUE, logfc.threshold = 0.25, test.use = 'bimod', min.pct = 0.1, random.seed =1)
View(DE_genes_cluster_4vsall)

# Save markers to file
write.csv(DE_genes_cluster_4vsall, 
          file = "~/Documents/Rachel's scRNAseq/48hr_DE_genes_cluster_4vsall.csv", 
          quote = FALSE, 
          row.names = TRUE)



#find conserved markers
Conserved_genes_2vs3 <- FindConservedMarkers(gcdata, ident.1 = 2, ident.2 = 3, grouping.var = "timepoint", only.pos = TRUE, min.diff.pct = 0.25, min.pct = 0.25, logfc.threshold = 0.25)
View(Conserved_genes_2vs3)

library(tidyverse)
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')



#find all markers for all clusters
combined_markers <- FindAllMarkers(object = gcdata, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25)  
View(combined_markers) 

save(day2, file = "~/Documents/Rachel_scRNAseq/gcdata_48hr_final_20220830.RData")


#label transfer method to compare our cells to fetal brain reference dataset
library(Seurat)
library(ggplot2)

pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = stri_11pcw, query = gcdata,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = AdultNeuroExNac_V2$EmbryoAdultNuclei,
                            dims = 1:30, k=5)
gcdata <- AddMetaData(gcdata, metadata = predictions)
gcdata$prediction.match <- gcdata$predicted.id == gcdata$SCINA_celltypes_3
  table(gcdata$prediction.match)
table(pancreas.query$predicted.id)


#label transfer method 2

gcdata <- RunPCA(gcdata)
gcdata <- RunUMAP(gcdata, assay = "RNA", dims = 1:20)
gcdata <- FindNeighbors(gcdata, dims = 1:15)
gcdata <- FindClusters(gcdata, dims = 1:15, resolution = 0.5)
gcdata<- RunUMAP(gcdata, assay = "RNA", dims = 1:15)
DimPlot(gcdata, reduction = "umap", group.by = "seurat_clusters", pt.size = .25) 
g <- FindAllMarkers(gcdata, assay = "RNA")


avg_expr_ref <- sapply(sort(unique(stri_11pcw$Cell_type)), function(ct) rowMeans(stri_11pcw@assays$RNA@data[,which(stri_11pcw$Cell_type == ct)] ))
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


#use gene module score to calculate similarity to striatum
library(RColorBrewer)
stri <- read.csv("~/Desktop/striatum_DE_genes.csv")
stri <- as.matrix(stri)

gcdata <- AddModuleScore(gcdata,
                       features = list(stri),
                       name="striatum")

# Plot scores
plot2 <- FeaturePlot(gcdata,
            features = "all_relevant_stri_celltypes1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

plot1 <- UMAPPlot(gcdata, group.by="seurat_clusters", label=T)
plot1 | plot2

FeaturePlot(gcdata,
            features = "GAD.interneuron.2", label = TRUE, repel = TRUE, cols = c("white", "red")) 

save(gcdata, file = "~/Documents/week4_final_w_genemodulescores.RData")


#convert h5ad to seurat
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

Convert("~/Downloads/Bocchi_Science2021_Striatum.h5ad", dest = "h5seurat", overwrite = TRUE)
striatum <- LoadH5Seurat("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/Bocchi_Science2021_Striatum.h5seurat", assays = "RNA")
stri_11pcw <- subset(x = striatum, subset = (PCW_._DAY == "10PCW6D" | PCW_._DAY == "11PCW3D"))

stri_11pcw_slc32a1_only <- subset(stri_11pcw, subset = SLC32A1 >0)
counts_matrix <- GetAssayData(stri_11pcw_slc32a1_only, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0("~/Desktop/", 'stri_11pcw_slc32a1_only_counts.mtx'))
df <- as.data.frame(counts_matrix)
write.csv(df, file='~/Desktop/stri_11pcw_slc32a1_only_counts.csv', quote=F, row.names=TRUE)

stri_11PCW3D_slc32a1 <- subset(stri_11pcw_slc32a1_only, subset = PCW_._DAY == "11PCW3D")
stri_10PCW6D_slc32a1 <- subset(stri_11pcw_slc32a1_only, subset = PCW_._DAY == "10PCW6D")

counts_matrix <- GetAssayData(stri_10PCW6D_slc32a1, assay='RNA', slot='counts')
df <- as.data.frame(counts_matrix)
write.csv(df, file='~/Desktop/stri_10PCW6D_slc32a1_only_counts.csv', quote=F, row.names=TRUE)




  
  #label transfer method for all GW14 fetal brain regions 

  load(file = "~/Desktop/GW14_allregions_merged_20231113.RData")
  
  avg_expr_ref <- sapply(sort(unique(combined$region.ident)), function(ct) rowMeans(combined@assays$RNA@data[,which(combined$region.ident == ct)] ))
  
  avg_expr_ds1 <- sapply(levels(gcdata@active.ident), function(ct) rowMeans(gcdata@assays$RNA@data[,which(gcdata@active.ident == ct)]))
  genes2cor <- intersect(VariableFeatures(combined), rownames(gcdata))
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
  
  anchors <- FindTransferAnchors(reference = combined, query = gcdata, dims = 1:30, npcs = 30)
  predictions <- TransferData(anchorset = anchors, refdata = combined$region.ident, dims = 1:30, k.weight = 3)
  gcdata$celltype_transfer <- predictions$predicted.id
  
  plot1 <- UMAPPlot(gcdata, label=T)
  plot2 <- UMAPPlot(gcdata, group.by="celltype_transfer", label=T)
  plot1 | plot2
  
  UMAPPlot(All.integrated, group.by="sampletype", label=T)
  
  
  
  #integration of striatum and hypo datasets
  
  # split the merged  dataset into a list  seurat objects 
  ifnb.list <- SplitObject(combined, split.by = "region.ident")
  
  
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
  GW20_brain_regions_combined <- IntegrateData(anchorset = immune.anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(GW20_brain_regions_integrated_slc32a1_only) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  GW20_brain_regions_integrated_slc32a1_only <- ScaleData(GW20_brain_regions_integrated_slc32a1_only, verbose = FALSE)
  GW20_brain_regions_integrated_slc32a1_only <- FindVariableFeatures(GW20_brain_regions_integrated_slc32a1_only, selection.method = "vst", nfeatures = 2000)
  GW20_brain_regions_integrated_slc32a1_only <- RunPCA(GW20_brain_regions_integrated_slc32a1_only, npcs = 30, verbose = FALSE)
  GW20_brain_regions_integrated_slc32a1_only <- RunUMAP(GW20_brain_regions_integrated_slc32a1_only, reduction = "pca", dims = 1:30)
  GW20_brain_regions_integrated_slc32a1_only <- FindNeighbors(GW20_brain_regions_integrated_slc32a1_only, reduction = "pca", dims = 1:30)
  GW20_brain_regions_integrated_slc32a1_only <- FindClusters(GW20_brain_regions_integrated_slc32a1_only, resolution = 0.5)
  
  # Visualization
  p1 <- DimPlot(GW20_brain_regions_integrated_slc32a1_only, reduction = "umap", group.by = "region.ident")
  p2 <- DimPlot(GW20_brain_regions_integrated_slc32a1_only, reduction = "umap", label = TRUE, repel = TRUE)
  p1 + p2
  
  
  
  # To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
  DimPlot(GW20_hypo_stri_regions_combined, reduction = "umap", split.by = "stim")
  
  DefaultAssay(GW20_brain_regions_integrated_slc32a1_only) <- "RNA"
  
  
  #pull of interneruons
  GW20_brain_regions_integrated_slc32a1_only <- subset(GW20_brain_regions_combined, subset = SLC32A1 >0)
  save(GW20_brain_regions_integrated_slc32a1_only, file = "~/Desktop/GW20_brain_regions_integrated_slc32a1_only.RData")
  gcdata_slc32a1_only <- subset(gcdata, subset = SLC32A1 >0)
  save(gcdata_slc32a1_only, file = "~/Desktop/week4_slc32a1_only.RData")
  
  
  DimPlot(GW20_hypo_stri_regions_combined_slc32a1_only, reduction = "umap", group.by = "region.ident")
  
  GW20_brain_regions_integrated_slc32a1_only <- FindVariableFeatures(GW20_brain_regions_integrated_slc32a1_only, selection.method = "vst", nfeatures = 2000)
  
  #label transfer method
  
  avg_expr_ref <- sapply(sort(unique(GW20_brain_regions_integrated_slc32a1_only$region.ident)), function(ct) rowMeans(GW20_brain_regions_integrated_slc32a1_only@assays$RNA@data[,which(GW20_brain_regions_integrated_slc32a1_only$region.ident == ct)] ))
  
  avg_expr_ds1 <- sapply(levels(gcdata_slc32a1_only@active.ident), function(ct) rowMeans(gcdata_slc32a1_only@assays$RNA@data[,which(gcdata_slc32a1_only@active.ident == ct)]))
  genes2cor <- intersect(VariableFeatures(GW20_brain_regions_integrated_slc32a1_only), rownames(gcdata_slc32a1_only))
  genes2cor <- intersect(rownames(GW20_brain_regions_integrated_slc32a1_only), rownames(gcdata_slc32a1_only))
  corr2ref_cl <- cor(avg_expr_ds1[genes2cor,], avg_expr_ref[genes2cor,], method="pearson")
  
  library(gplots)
  heatmap.2(corr2ref_cl, scale="none", trace="none", dendrogram = "none", key=T, keysize=1, key.par=list(mar=c(3,3.5,2,0)), key.title = 'Correlation',
            key.xlab = 'Value', margins=c(14,16),
            labRow = colnames(avg_expr_ds1), labCol = colnames(avg_expr_ref), cexRow=0.8, cexCol=0.8,
            col=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")))(30))
  legend("right", title = "Tissues",legend=c("group1","group2"), 
         fill=c("red","green"), cex=0.8, box.lty=0)
  
  #transcriptome similarity on cellular level
  ranked_expr_ref <- apply(avg_expr_ref[genes2cor,],2,rank)
  library(presto)
  ranked_expr_ds1 <- rank_matrix(gcdata_slc32a1_only@assays$RNA@data[genes2cor,])$X_ranked
  library(qlcMatrix)
  corr2ref_cell <- corSparse(ranked_expr_ds1, ranked_expr_ref)
  ct_maxcor <- colnames(avg_expr_ref)[apply(corr2ref_cell, 1, which.max)]
  gcdata_slc32a1_only$celltype_maxcor <- ct_maxcor
  plot1 <- UMAPPlot(gcdata_slc32a1_only, label=T)
  plot2 <- UMAPPlot(res, group.by="type", label=F)
  plot1 | plot2

  
  save(combined, file = "~/Desktop/GW20_caudate_NA_putamen_hypo_merged.RData")
  save(GW20_hypo_stri_regions_combined, file = "~/Desktop/GW20_caudate_NA_putamen_hypo_inegrated_badhuri.RData")
  
  
  #identify DE genes across conditions
  
  #change identity
  Idents(object = GW20_hypo_stri_regions_combined_slc32a1_only) <- GW20_hypo_stri_regions_combined_slc32a1_only@meta.data$region.ident
  caudate_genes <- FindAllMarkers(GW20_hypo_stri_regions_combined_slc32a1_only)
  View(caudate_genes)
  
  write.csv(caudate_genes, file = "~/Desktop/hypo_DE_genes_badhuri_findallmarkers_top10.txt")

  hypo <- read.delim("~/Desktop/findallmarkers_putamen_DE_genes_GW20_badhuri_top500.txt", sep = "\t", stringsAsFactors = FALSE)
  
  
  gcdata_slc32a1_only <- AddModuleScore(gcdata_slc32a1_only,
                         features = hypo,
                         name="putamen_top500")
  
  # Plot scores
  plot2 <- FeaturePlot(gcdata_slc32a1_only,
                       features = "hypo_top2001", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  plot3 <- FeaturePlot(gcdata_slc32a1_only,
                       features = "caudate1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  plot4 <- FeaturePlot(gcdata_slc32a1_only,
                       features = "putamen_top5001", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  plot5 <- FeaturePlot(gcdata_slc32a1_only,
                       features = "hypothalamus1", label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  plot1 <- UMAPPlot(gcdata_slc32a1_only, group.by="seurat_clusters", label=T)
  par(mfrow = c(3, 3))
  
  plot1 + plot2 + plot3 + plot4 + plot5

  plot1 + plot4
  
  
  week4_interneurons_only <- subset(gcdata, SLC32A1 > 0 & GAD1 > 0 & GAD2 > 0)
  t <- subset(gcdata, subset = GABA_mod.score1 >0)
  gc_48hr <- subset(gcdata, subset = timepoint == "48h_pD")
  
  features <- list(c("SLC32A1", "GAD1", "GAD2"))
  gcdata <- AddModuleScore(gcdata,
                                        features = "SLC32A1",
                                        name="GABA_mod.score")
  
FeaturePlot(gcdata,
                       features = "GABA_mod.score1", label = TRUE, repel = TRUE)

sum(GetAssayData(object = gcdata, slot = "counts")["GAD1",] >0 )

length(WhichCells(object = gcdata, expression = SLC32A1 > 0 & GAD1 > 0 & GAD2))/nrow(gcdata@meta.data)*100
length(WhichCells(object = gcdata, slot = "data", expression = SLC32A1 > 0 & GAD1 > 0 & GAD2 > 0))


devtools::install_github("prestodb/RPresto")
library(RPresto)

g <- FindAllMarkers(gcdata, logfc.threshold = 0.25, slot = "scale.data")

g %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_diff) -> top10
DoHeatmap(gcdata, features = top10$gene)

load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_LB.RData")


#integration of all AD scRNA-seq timepoints using seurat

gcdata <- NormalizeData(gcdata)
gcdata <- FindVariableFeatures(gcdata, selection.method = "vst", nfeatures = 2000)

gcdata <- ScaleData(gcdata, verbose = FALSE)
gcdata <- RunPCA(gcdata, npcs = 30, verbose = FALSE)
gcdata <- RunUMAP(gcdata, reduction = "pca", dims = 1:30)
gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:30)
gcdata <- FindClusters(gcdata, resolution = 0.5)

# split the merged  dataset into a list  seurat objects 
ifnb.list <- SplitObject(gcdata, split.by = "timepoint")


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list, dims = 1:30)

#perform integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
gcdata_seurat_integrated <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(gcdata_seurat_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
gcdata_seurat_integrated <- ScaleData(gcdata_seurat_integrated, verbose = FALSE)
gcdata_seurat_integrated <- FindVariableFeatures(gcdata_seurat_integrated, selection.method = "vst", nfeatures = 2000)
gcdata_seurat_integrated <- RunPCA(gcdata_seurat_integrated, npcs = 30, verbose = FALSE)
gcdata_seurat_integrated <- RunUMAP(gcdata_seurat_integrated, reduction = "pca", dims = 1:30)
gcdata_seurat_integrated <- RunUMAP(gcdata_seurat_integrated, dims = 1:30)
gcdata_seurat_integrated <- FindNeighbors(gcdata_seurat_integrated)
gcdata_seurat_integrated <- FindClusters(gcdata_seurat_integrated, resolution = 0.6)

# Visualization
p1 <- DimPlot(gcdata_seurat_integrated, reduction = "umap", group.by = "timepoint")
p2 <- DimPlot(gcdata_seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

save(gcdata_seurat_integrated, file = "~/Desktop/gcdata_AD_seurat_integrated.RData")



#handling of individual timepoints

load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_48hr_final_20220830.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_day7_final_20220829.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_week4_celltype_annotations_final.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_1wk_2wk_3wk.RData")
load("~/Desktop/wk1_gcdata.RData")
load("~/Desktop/wk2_gcdata.RData")
load("~/Desktop/wk3_gcdata.RData")
load("~/Desktop/gcdata_all_timepoints_CC_genes_regressed.RData")

wk1 <- subset(gc_1wk_2wk_3wk, subset = timepoint == "1wk_glia")
wk2 <- subset(gc_1wk_2wk_3wk, subset = timepoint == "2wk_glia")
wk3 <- subset(gc_1wk_2wk_3wk, subset = timepoint == "3wk_glia")
save(wk1, file = "~/Desktop/wk1_gcdata.RData")
save(wk2, file = "~/Desktop/wk2_gcdata.RData")
save(wk3, file = "~/Desktop/wk3_gcdata.RData")

load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_LB.RData")

wk1 <- subset(gcdata_mito_filtered, subset = timepoint == "1wk_glia")
wk2 <- subset(gcdata_mito_filtered, subset = timepoint == "2wk_glia")
wk3 <- subset(gcdata_mito_filtered, subset = timepoint == "3wk_glia")
wk4 <- subset(gcdata_mito_filtered, subset = timepoint == "4wk_glia")
day7 <- subset(gcdata_mito_filtered, subset = timepoint == "7d_postDox")
day2 <- subset(gcdata_mito_filtered, subset = timepoint == "48h_pD")
dualSMAD <- subset(gcdata_mito_filtered, subset = timepoint == "dualSMAD")

wk1 <- FindVariableFeatures(wk1, selection.method = "vst", nfeatures = 1000)
wk2 <- FindVariableFeatures(wk2, selection.method = "vst", nfeatures = 1000)
wk3 <- FindVariableFeatures(wk3, selection.method = "vst", nfeatures = 1000)
wk4 <- FindVariableFeatures(wk4, selection.method = "vst", nfeatures = 1000)
day7 <- FindVariableFeatures(day7, selection.method = "vst", nfeatures = 1000)
day2 <- FindVariableFeatures(day2, selection.method = "vst", nfeatures = 1000)
dualSMAD <- FindVariableFeatures(dualSMAD, selection.method = "vst", nfeatures = 1000)
gcdata_mito_filtered <- FindVariableFeatures(gcdata_mito_filtered, selection.method = "vst", nfeatures = 1000)

write.csv(day2@assays$RNA@var.features, file = "~/Desktop/day2_varfeatures_1000.csv")
write.csv(day7@assays$RNA@var.features, file = "~/Desktop/day7_varfeatures_1000.csv")
write.csv(wk4@assays$RNA@var.features, file = "~/Desktop/week4_varfeatures_1000.csv")
write.csv(gcdata_mito_filtered@assays$RNA@var.features, file = "~/Desktop/alltimepoints_varfeatures_1000.csv")
write.csv(wk1@assays$RNA@var.features, file = "~/Desktop/week1_varfeatures_1000.csv")
write.csv(wk2@assays$RNA@var.features, file = "~/Desktop/week2_varfeatures_1000.csv")
write.csv(wk3@assays$RNA@var.features, file = "~/Desktop/week3_varfeatures_1000.csv")
write.csv(dualSMAD@assays$RNA@var.features, file = "~/Desktop/dualSMAD_varfeatures_1000.csv")


#cell cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
gcdata_mito_filtered_all_10_wk2_15 <- NormalizeData(gcdata_mito_filtered_all_10_wk2_15)
gcdata_mito_filtered_all_10_wk2_15 <- FindVariableFeatures(gcdata_mito_filtered_all_10_wk2_15, selection.method = "vst", nfeatures = 2000)
gcdata_mito_filtered_all_10_wk2_15 <- ScaleData(gcdata_mito_filtered_all_10_wk2_15, features = rownames(gcdata_mito_filtered_all_10_wk2_15))

gcdata_mito_filtered_all_10_wk2_15 <- CellCycleScoring(gcdata_mito_filtered_all_10_wk2_15, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(gcdata[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(gcdata_mito_filtered_all_10_wk2_15, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "timepoint")

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
gcdata_mito_filtered_all_10_wk2_15 <- RunPCA(gcdata_mito_filtered_all_10_wk2_15, features = c(s.genes, g2m.genes))
DimPlot(gcdata_mito_filtered_all_10_wk2_15, reduction = "pca")

#Regress out cell cycle scores during data scaling
gcdata_mito_filtered_all_10_wk2_15_CC_regressed <- ScaleData(gcdata_mito_filtered_all_10_wk2_15, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(gcdata_mito_filtered_all_10_wk2_15))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
gcdata_mito_filtered_all_10_wk2_15_CC_regressed <- RunPCA(gcdata_mito_filtered_all_10_wk2_15_CC_regressed, features = VariableFeatures(gcdata_mito_filtered_all_10_wk2_15_CC_regressed), nfeatures.print = 10)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
gcdata_mito_filtered_all_10_wk2_15_CC_regressed <- RunPCA(gcdata_mito_filtered_all_10_wk2_15_CC_regressed, features = c(s.genes, g2m.genes))
DimPlot(gcdata_mito_filtered_all_10_wk2_15_CC_regressed)
DimPlot(gcdata_mito_filtered_all_10_wk2_15_CC_regressed, reduction = "pca", group.by = "Phase")

save(gcdata_mito_filtered_all_10_wk2_15_CC_regressed, file = "~/Desktop/gcdata_mito_10_wk2_15_filtered_CC_regressed_20230817.RData")

counts_matrix <- GetAssayData(gcdata_CC_regressed_CSS_integrated, assay='RNA', slot='counts')
df <- as.matrix(counts_matrix)
write.csv(df, "~/Desktop/CSS_integratedCC_regressed_RNA_counts.csv")
TFs <-c("ARID3B", "ARID3C", "ZBED6", "AHRR", "ARNTL2", "ASCL3", "ASCL4")

TF1 <- FetchData(gcdata_CC_regressed_CSS_integrated, vars = TFs1)
TF2 <- FetchData(gcdata_CC_regressed_CSS_integrated, vars = TFs)
TF3 <- FetchData(gcdata_CC_regressed_CSS_integrated, vars = Tfs2)
TF4 <- FetchData(gcdata_CC_regressed_CSS_integrated, vars = Tfs3)
TF5 <- FetchData(gcdata_CC_regressed_CSS_integrated, vars = Tfs4)
TF6 <- FetchData(gcdata_CC_regressed_CSS_integrated, vars = Tfs5)
write.csv(TF1, "~/Desktop/TF1_RNA_data.csv")
write.csv(TF2, "~/Desktop/TF2_RNA_data.csv")
write.csv(TF3, "~/Desktop/TF3_RNA_data.csv")
write.csv(TF4, "~/Desktop/TF4_RNA_data.csv")
write.csv(TF5, "~/Desktop/TF5_RNA_data.csv")
write.csv(TF6, "~/Desktop/TF6_RNA_data.csv")
TFs <- read.csv("~/Desktop/All_TFs_RNA_data.csv")

gabascore = gcdata_CC_regressed_CSS_integrated$GAD.interneuron
df <- as.matrix(gabascore)
write.csv(df, "~/Desktop/Gabascores.csv")
TFs <- read.csv("~/Desktop/TF_exp_sum_for_corr_gabascore.csv", row.names = 1)
TF <- read.csv("~/Desktop/All_TFs_RNA_data.csv")
cor(TFs, df, method = "pearson")

matrix<-gcdata_CC_regressed_CSS_integrated@assays$RNA@counts
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(gcdata_CC_regressed_CSS_integrated@meta.data$GAD.interneuron)
correlations<-apply(matrix_mod,1,function(x){cor(gene,x, method = "pearson")})
write.csv(correlations, "~/Desktop/all_genes_correlation_gabascore_rna_counts.csv")
gcdata_CC_regressed_CSS_integrated
TFs1 <-c("ARID3B", "ARID3C", "ZBED6", "AHRR", "ARNTL2", "ASCL3", "ASCL4", "ASCL5", "ATOH8", "HELT", "HES3", "HES4", "HEYL", "HIF3A", "LYL1", "MLXIP", "MXD1", "MXD3", "MXD4", "MYCL", "NEUROD4", "NEUROD6", "NEUROG3", "NPAS1", "NPAS3", "NPAS4", "SCX", "SIM1", "SIM2", "TAL2", "TCF15", "TCF23", "TCF24", "CREBL2", "CREBZF", "BNC1", "BNC2", "DPF3", "FEZF2", "HIVEP3", "IKZF4", "IKZF5", "MYT1L", "PRDM16", "SALL1", "SALL4", "SP6", "ST18", "ZBTB34", "ZFP2", "ZFP30", "ZFY", "ZNF117", "ZNF652", "ZNF664", "ZNF679", "ZNF691", "ZNF697", "ZNF705A", "ZNF705B", "ZNF705D", "ZNF852", "ZSCAN10", "ZEB2", "ZFHX4", "TRERF1", "JRKL", "TIGD2", "TIGD6", "LIN28A", "YBX2", "RBPJL", "SATB2", "TFDP2", "TFDP3", "FOXB2", "FOXD4", "FOXD4L1", "FOXD4L3", "FOXD4L5", "FOXD4L6", "FOXE3", "FOXI2", "FOXI3", "FOXK2", "FOXP4", "FOXS1", "HMG20B", "DBX1", "DBX2", "HHEX", "IRX6", "NKX2-6", "RHOXF2B", "SEBOX", "NAIF1", "TERF2", "NR1H3", "PROX2", "RFX6", "DEAF1", "TBX10", "TBPL2", "LIN54", "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D", "TFAP2E", "ARID3A", "ARID5A", "ARID5B", "HMGA1", "HMGA2", "ZBED1", "AHR", "ARNT", "ARNT2", "ARNTL", "ASCL1", "ASCL2", "ATOH1", "ATOH7", "BHLHA15", "BHLHE22", "BHLHE23", "BHLHE40", "BHLHE41", "CCDC169-SOHLH2", "CLOCK", "EPAS1", "FERD3L", "FIGLA", "HAND1", "HAND2", "HES1", "HES2", "HES5", "HES6", "HES7", "HEY1", "HEY2", "HIF1A", "MAX", "MESP1", "MESP2", "MITF", "MLX", "MLXIPL", "MNT", "MSC", "MSGN1", "MXI1", "MYC", "MYCN", "MYF5", "MYF6", "MYOD1", "MYOG", "NEUROD1", "NEUROD2", "NEUROG1", "NEUROG2", "NHLH1", "NHLH2", "NPAS2", "OLIG1", "OLIG2", "OLIG3", "PTF1A", "SOHLH2", "SREBF1", "SREBF2", "TAL1", "TCF12", "TCF21", "TCF3", "TCF4", "TCFL5", "TFAP4", "TFE3", "TFEB", "TFEC", "TWIST1", "TWIST2", "USF1", "USF2", "AC023509.3", "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6", "ATF6B", "ATF7", "BACH1", "BACH2", "BATF", "BATF3", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CREB1", "CREB3", "CREB3L1", "CREB3L2", "CREB3L4", "CREB5", "CREM", "DBP", "DDIT3", "FOS", "FOSB", "FOSL1", "FOSL2", "HLF", "JDP2", "JUN", "JUNB", "JUND", "MAF", "MAFA", "MAFB", "MAFF", "MAFG", "MAFK", "NFE2", "NFE2L1", "NFE2L2", "NFE2L3", "NFIL3", "NRL", "TEF", "XBP1", "AC138696.1", "BCL11A", "BCL11B", "BCL6", "BCL6B", "CTCF", "CTCFL", "DPF1", "E4F1", "EGR1", "EGR2", "EGR3", "EGR4", "FEZF1", "GFI1", "GFI1B", "GLI1", "GLI2", "GLI3", "GLI4", "GLIS1", "GLIS2", "GLIS3", "GTF3A", "GZF1", "HIC1", "HIC2", "HINFP", "HIVEP1", "HIVEP2", "HKR1", "IKZF1", "IKZF2", "IKZF3", "INSM1", "KLF1", "KLF10", "KLF11", "KLF12", "KLF13", "KLF14", "KLF15", "KLF16", "KLF17", "KLF2", "KLF3", "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9", "MAZ", "MECOM", "MTF1", "MYNN", "MZF1", "OSR1", "OSR2", "OVOL1", "OVOL2", "PLAG1", "PLAGL1", "PLAGL2", "PRDM1", "PRDM12", "PRDM14", "PRDM15", "PRDM4", "PRDM6", "PRDM9", "RBAK", "REST", "RREB1", "SALL2", "SCRT1", "SCRT2", "SNAI1", "SNAI2", "SNAI3", "SP1", "SP2", "SP3", "SP4", "SP5", "SP7", "SP8", "SP9", "VEZF1", "WT1", "YY1", "YY2", "ZBTB1", "ZBTB11", "ZBTB12", "ZBTB14", "ZBTB16", "ZBTB18", "ZBTB2", "ZBTB20", "ZBTB22", "ZBTB26", "ZBTB3", "ZBTB32", "ZBTB33", "ZBTB37", "ZBTB4", "ZBTB42", "ZBTB43", "ZBTB44", "ZBTB45", "ZBTB48", "ZBTB49", "ZBTB6", "ZBTB7A", "ZBTB7B", "ZBTB7C", "ZFP1", "ZFP14", "ZFP28", "ZFP3", "ZFP41", "ZFP42", "ZFP57", "ZFP64", "ZFP69", "ZFP69B", "ZFP82", "ZFP90", "ZFP92", "ZFX", "ZIC1", "ZIC2", "ZIC3", "ZIC4", "ZIC5", "ZIK1", "ZIM2", "ZIM3", "ZKSCAN1", "ZKSCAN2", "ZKSCAN3", "ZKSCAN5", "ZKSCAN7", "ZNF10", "ZNF100", "ZNF101", "ZNF114", "ZNF12", "ZNF121", "ZNF124", "ZNF132", "ZNF133", "ZNF134", "ZNF135", "ZNF136", "ZNF140", "ZNF141", "ZNF143", "ZNF146", "ZNF148", "ZNF154", "ZNF157", "ZNF16", "ZNF169", "ZNF17", "ZNF174", "ZNF175", "ZNF177", "ZNF18", "ZNF180", "ZNF181", "ZNF182", "ZNF184", "ZNF189", "ZNF19", "ZNF197", "ZNF2", "ZNF200", "ZNF202", "ZNF205", "ZNF211", "ZNF212", "ZNF213", "ZNF214", "ZNF217", "ZNF219", "ZNF22", "ZNF222", "ZNF223", "ZNF224", "ZNF225", "ZNF23", "ZNF232", "ZNF235", "ZNF236", "ZNF248", "ZNF25", "ZNF250", "ZNF254", "ZNF257", "ZNF26", "ZNF260", "ZNF263", "ZNF264", "ZNF266", "ZNF267", "ZNF273", "ZNF274", "ZNF276", "ZNF28", "ZNF280A", "ZNF281", "ZNF282", "ZNF283", "ZNF284", "ZNF285", "ZNF287", "ZNF296", "ZNF3", "ZNF30", "ZNF300", "ZNF302", "ZNF304", "ZNF311", "ZNF316", "ZNF317", "ZNF32", "ZNF320", "ZNF322", "ZNF324", "ZNF324B", "ZNF329", "ZNF331", "ZNF333", "ZNF334", "ZNF337", "ZNF33A", "ZNF33B", "ZNF34", "ZNF341", "ZNF343", "ZNF345", "ZNF35", "ZNF350", "ZNF354A", "ZNF354B", "ZNF354C", "ZNF37A", "ZNF382", "ZNF383", "ZNF384", "ZNF385D", "ZNF394", "ZNF396", "ZNF398", "ZNF41", "ZNF410", "ZNF415", "ZNF416", "ZNF417", "ZNF418", "ZNF419", "ZNF423", "ZNF425", "ZNF429", "ZNF430", "ZNF431", "ZNF432", "ZNF433", "ZNF436", "ZNF439", "ZNF44", "ZNF440", "ZNF441", "ZNF442", "ZNF443", "ZNF444", "ZNF445", "ZNF449", "ZNF45", "ZNF454", "ZNF460", "ZNF467", "ZNF468", "ZNF479", "ZNF480", "ZNF483", "ZNF484", "ZNF485", "ZNF486", "ZNF487", "ZNF490", "ZNF492", "ZNF496", "ZNF501", "ZNF502", "ZNF506", "ZNF513", "ZNF519", "ZNF521", "ZNF525", "ZNF527", "ZNF528", "ZNF529", "ZNF530", "ZNF534", "ZNF540", "ZNF543", "ZNF547", "ZNF548", "ZNF549", "ZNF550", "ZNF552", "ZNF554", "ZNF555", "ZNF557", "ZNF558", "ZNF561", "ZNF562", "ZNF563", "ZNF564", "ZNF565", "ZNF566", "ZNF567", "ZNF570", "ZNF571", "ZNF573", "ZNF574", "ZNF580", "ZNF581", "ZNF582", "ZNF584", "ZNF585A", "ZNF586", "ZNF587", "ZNF589", "ZNF594", "ZNF595", "ZNF596", "ZNF597", "ZNF605", "ZNF610", "ZNF611", "ZNF613", "ZNF614", "ZNF615", "ZNF616", "ZNF619", "ZNF620", "ZNF621", "ZNF626", "ZNF627", "ZNF628", "ZNF641", "ZNF649", "ZNF655", "ZNF658", "ZNF660", "ZNF662", "ZNF667", "ZNF669", "ZNF671", "ZNF674", "ZNF675", "ZNF677", "ZNF680", "ZNF681", "ZNF682", "ZNF684", "ZNF69", "ZNF692", "ZNF695", "ZNF7", "ZNF701", "ZNF704", "ZNF705G", "ZNF707", "ZNF708", "ZNF71", "ZNF711", "ZNF713", "ZNF714", "ZNF716", "ZNF718", "ZNF727", "ZNF730", "ZNF735", "ZNF736", "ZNF737", "ZNF74", "ZNF740", "ZNF749", "ZNF75A", "ZNF75D", "ZNF76", "ZNF761", "ZNF764", "ZNF765", "ZNF766", "ZNF768", "ZNF77", "ZNF770", "ZNF771", "ZNF774", "ZNF776", "ZNF777", "ZNF778", "ZNF780A", "ZNF782", "ZNF783", "ZNF784", "ZNF785", "ZNF786", "ZNF787", "ZNF789", "ZNF79", "ZNF790", "ZNF791", "ZNF792", "ZNF793", "ZNF799", "ZNF8", "ZNF805", "ZNF808", "ZNF81", "ZNF816", "ZNF821", "ZNF823", "ZNF84", "ZNF846", "ZNF85", "ZNF860", "ZNF879", "ZNF880", "ZNF883", "ZNF891", "ZNF90", "ZNF93", "ZNF98", "ZSCAN1", "ZSCAN16", "ZSCAN22", "ZSCAN23", "ZSCAN26", "ZSCAN29", "ZSCAN30", "ZSCAN31", "ZSCAN4", "ZSCAN5A", "ZSCAN5C", "ZSCAN9", "PATZ1", "ZNF524", "ZNF653", "ZEB1", "ZFHX3", "ZNF541", "NFYA", "CENPB", "CENPBD1", "TIGD1", "LIN28B", "YBX1", "YBX3", "RBPJ", "CUX1", "CUX2", "ONECUT1", "ONECUT2", "ONECUT3", "SATB1", "CXXC1", "CXXC5", "DNMT1", "KDM2B", "TET1", "KMT2A", "DMRT1", "DMRT2", "DMRT3", "DMRTA1", "DMRTA2", "DMRTC2", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8", "TFDP1", "EBF1", "EBF2", "EBF3", "EBF4", "EHF", "ELF1", "ELF2", "ELF4", "ELF5", "ELK1", "ELK3", "ELK4", "ERF", "ERG", "ETS1", "ETS2", "ETV1", "ETV2", "ETV3", "ETV3L", "ETV4", "ETV5", "ETV6", "ETV7", "FEV", "FLI1", "GABPA", "SPDEF", "SPI1", "SPIB", "SPIC", "ELF3", "FOXA1", "FOXA2", "FOXA3", "FOXB1", "FOXC1", "FOXC2", "FOXD1", "FOXD2", "FOXD3", "FOXD4L4", "FOXE1", "FOXF1", "FOXF2", "FOXG1", "FOXH1", "FOXI1", "FOXJ1", "FOXJ2", "FOXJ3", "FOXK1", "FOXL1", "FOXL2", "FOXM1", "FOXN1", "FOXN2", "FOXN3", "FOXN4", "FOXO1", "FOXO3", "FOXO4", "FOXO6", "FOXP1", "FOXP2", "FOXP3", "FOXQ1", "FOXR1", "FOXR2", "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6", "TRPS1", "GCM1", "GCM2", "GRHL1", "GRHL2", "TFCP2", "TFCP2L1", "UBP1", "GTF2I", "GTF2IRD1", "BBX", "CIC", "HBP1", "LEF1", "SOX1", "SOX10", "SOX11", "SOX12", "SOX13", "SOX14", "SOX15", "SOX17", "SOX18", "SOX2", "SOX21", "SOX3", "SOX30", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8", "SOX9", "SRY", "TCF7", "TCF7L1", "TCF7L2", "ALX1", "ALX3", "ALX4", "ANHX", "ARGFX", "ARX", "BARHL1", "BARHL2", "BARX1", "BARX2", "BSX", "CDX1", "CDX2", "CDX4", "CRX", "DLX1", "DLX2", "DLX3", "DLX4", "DLX5", "DLX6", "DMBX1", "DPRX", "DRGX", "DUX1", "DUX3", "DUX4", "DUXA", "EMX1", "EMX2", "EN1", "EN2", "ESX1", "EVX1", "EVX2", "GBX1", "GBX2", "GSC", "GSC2", "GSX1", "GSX2", "HDX", "HESX1", "HLX", "HMBOX1", "HMX1", "HMX2", "HMX3", "HNF1A", "HNF1B", "HOMEZ", "HOXA1", "HOXA10", "HOXA11", "HOXA13", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXB1", "HOXB13", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9", "HOXC10", "HOXC11", "HOXC12", "HOXC13", "HOXC4", "HOXC5", "HOXC6", "HOXC8", "HOXC9", "HOXD1", "HOXD10", "HOXD11", "HOXD12", "HOXD13", "HOXD3", "HOXD4", "HOXD8", "HOXD9", "IRX1", "IRX2", "IRX3", "IRX4", "IRX5", "ISL1", "ISL2", "ISX", "LBX1", "LBX2", "LHX1", "LHX2", "LHX3", "LHX4", "LHX5", "LHX6", "LHX8", "LHX9", "LMX1A", "LMX1B", "MEIS1", "MEIS2", "MEIS3", "MEOX1", "MEOX2", "MIXL1", "MNX1", "MSX1", "MSX2", "NANOG", "NANOGP8", "NKX1-1", "NKX1-2", "NKX2-1", "NKX2-2", "NKX2-3", "NKX2-4", "NKX2-5", "NKX2-8", "NKX3-1", "NKX3-2", "NKX6-1", "NKX6-2", "NKX6-3", "NOBOX", "NOTO", "OTP", "OTX1", "OTX2", "PBX1", "PBX2", "PBX3", "PBX4", "PDX1", "PHOX2A", "PHOX2B", "PITX1", "PITX2", "PITX3", "PKNOX1", "PKNOX2", "PROP1", "PRRX1", "PRRX2", "RAX", "RAX2", "RHOXF1", "RHOXF2", "SHOX", "SHOX2", "SIX1", "SIX2", "SIX3", "SIX4", "SIX5", "SIX6", "TGIF1", "TGIF2", "TGIF2LX", "TGIF2LY", "TLX1", "TLX2", "TLX3", "UNCX", "VAX1", "VAX2", "VENTX", "VSX1", "VSX2", "ZFHX2", "ZHX1", "PAX2", "PAX3", "PAX4", "PAX6", "PAX7", "POU1F1", "POU2F1", "POU2F2", "POU2F3", "POU3F1", "POU3F2", "POU3F3", "POU3F4", "POU4F1", "POU4F2", "POU4F3", "POU5F1", "POU5F1B", "POU6F1", "POU6F2", "HSF1", "HSF2", "HSF4", "HSF5", "HSFY1", "HSFY2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "MSANTD3", "BORCS8-MEF2B", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "SRF", "MBD2", "SETDB1", "MECP2", "MTERF1", "CDC5L", "MYB", "MYBL1", "MYBL2", "ZZZ3", "MYRF", "AR", "ESR1", "ESR2", "ESRRA", "ESRRB", "ESRRG", "HNF4A", "HNF4G", "NR0B1", "NR1D1", "NR1D2", "NR1H2", "NR1H4", "NR1I2", "NR1I3", "NR2C1", "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1", "PGR", "PPARA", "PPARD", "PPARG", "RARA", "RARB", "RARG", "RORA", "RORB", "RORC", "RXRA", "RXRB", "RXRG", "THRA", "THRB", "VDR", "BPTF", "CARF", "CEBPZ", "CPEB1", "GTF2B", "LTF", "MTF2", "NFE4", "NFYB", "NFYC", "NRF1", "PHF1", "POU2AF1", "PURA", "SKOR1", "SPZ1", "TOPORS", "XPA", "TP53", "TP63", "TP73", "PAX1", "PAX5", "PAX8", "PAX9", "LCOR", "LCORL", "PROX1", "NFAT5", "NFATC1", "NFATC2", "NFATC3", "NFATC4", "NFKB1", "NFKB2", "REL", "RELA", "RELB", "RFX1", "RFX2", "RFX3", "RFX4", "RFX5", "RFX7", "RUNX1", "RUNX2", "RUNX3", "AIRE", "GMEB1", "GMEB2", "SKOR2", "NFIA", "NFIB", "NFIC", "NFIX", "SMAD1", "SMAD3", "SMAD4", "SMAD5", "SMAD9", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "EOMES", "MGA", "T", "TBR1", "TBX1", "TBX15", "TBX18", "TBX19", "TBX2", "TBX20", "TBX21", "TBX22", "TBX3", "TBX4", "TBX5", "TBX6", "TBP", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "THAP1", "THAP12", "KDM5B", "ARID2", "AHCTF1", "AHDC1", "AKNA", "ASH1L", "CBX2", "DNTTIP1", "DOT1L", "GLYR1", "PHF20", "PHF21A", "PRR12", "SCML4", "SETBP1", "SRCAP", "C11orf95", "FAM200B", "SGSM2", "ZBED2", "ZBED3", "ZBED4", "ZBED5", "ZBED9", "BHLHA9", "NCOA1", "NCOA2", "NCOA3", "SOHLH1", "USF3", "POGK", "BATF2", "CREB3L3", "AC008770.3", "AC092835.1", "AEBP2", "AKAP8", "AKAP8L", "ANKZF1", "ATMIN", "CASZ1", "CCDC17", "CHAMP1", "CPXCR1", "DZIP1", "EEA1", "FAM170A", "FIZ1", "INSM2", "JAZF1", "KAT7", "KCMF1", "KIN", "L3MBTL1", "L3MBTL3", "L3MBTL4", "MYT1", "OVOL3", "PEG3", "PRDM10", "PRDM13", "PRDM2", "PRDM5", "PRDM8", "PRMT3", "RBSN", "REPIN1", "RLF", "SALL3", "SLC2A4RG", "TRAFD1", "TSHZ1", "TSHZ2", "TSHZ3", "WIZ", "ZBTB10", "ZBTB17", "ZBTB21", "ZBTB25", "ZBTB38", "ZBTB39", "ZBTB40", "ZBTB41", "ZBTB46", "ZBTB47", "ZBTB5", "ZBTB8A", "ZBTB8B", "ZBTB9", "ZFAT", "ZFP37", "ZFP62", "ZFP91", "ZFPM1", "ZFPM2", "ZKSCAN4", "ZKSCAN8", "ZMAT1", "ZMAT4", "ZNF107", "ZNF112", "ZNF131", "ZNF138", "ZNF14", "ZNF142", "ZNF155", "ZNF160", "ZNF165", "ZNF195", "ZNF20", "ZNF207", "ZNF208", "ZNF215", "ZNF221", "ZNF226", "ZNF227", "ZNF229", "ZNF230", "ZNF233", "ZNF234", "ZNF239", "ZNF24", "ZNF251", "ZNF253", "ZNF256", "ZNF268", "ZNF275", "ZNF280B", "ZNF280C", "ZNF280D", "ZNF286A", "ZNF286B", "ZNF292", "ZNF318", "ZNF319", "ZNF326", "ZNF335", "ZNF346", "ZNF347", "ZNF358", "ZNF362", "ZNF365", "ZNF366", "ZNF367", "ZNF385A", "ZNF385B", "ZNF385C", "ZNF391", "ZNF395", "ZNF397", "ZNF404", "ZNF407", "ZNF408", "ZNF414", "ZNF420", "ZNF426", "ZNF428", "ZNF43", "ZNF438", "ZNF446", "ZNF451", "ZNF461", "ZNF462", "ZNF469", "ZNF470", "ZNF471", "ZNF473", "ZNF474", "ZNF48", "ZNF488", "ZNF491", "ZNF493", "ZNF497", "ZNF500", "ZNF503", "ZNF507", "ZNF510", "ZNF511", "ZNF512B", "ZNF514", "ZNF516", "ZNF517", "ZNF518A", "ZNF518B", "ZNF526", "ZNF532", "ZNF536", "ZNF544", "ZNF546", "ZNF551", "ZNF556", "ZNF559", "ZNF560", "ZNF568", "ZNF569", "ZNF57", "ZNF572", "ZNF575", "ZNF576", "ZNF577", "ZNF578", "ZNF579", "ZNF583", "ZNF585B", "ZNF587B", "ZNF592", "ZNF598", "ZNF599", "ZNF600", "ZNF606", "ZNF607", "ZNF608", "ZNF609", "ZNF618", "ZNF623", "ZNF624", "ZNF625", "ZNF629", "ZNF630", "ZNF639", "ZNF644", "ZNF645", "ZNF646", "ZNF648", "ZNF654", "ZNF66", "ZNF665", "ZNF668", "ZNF670", "ZNF672", "ZNF676", "ZNF678", "ZNF683", "ZNF687", "ZNF688", "ZNF689", "ZNF696", "ZNF699", "ZNF70", "ZNF700", "ZNF703", "ZNF705E", "ZNF706", "ZNF709", "ZNF710", "ZNF717", "ZNF721", "ZNF724", "ZNF726", "ZNF728", "ZNF729", "ZNF732", "ZNF746", "ZNF747", "ZNF750", "ZNF763", "ZNF772", "ZNF773", "ZNF775", "ZNF780B", "ZNF781", "ZNF788", "ZNF80", "ZNF800", "ZNF804A", "ZNF804B", "ZNF813", "ZNF814", "ZNF827", "ZNF829", "ZNF83", "ZNF830", "ZNF831", "ZNF835", "ZNF836", "ZNF837", "ZNF841", "ZNF843", "ZNF844", "ZNF845", "ZNF850", "ZNF853", "ZNF865", "ZNF878", "ZNF888", "ZNF91", "ZNF92", "ZNF99", "ZSCAN12", "ZSCAN18", "ZSCAN2", "ZSCAN20", "ZSCAN21", "ZSCAN25", "ZSCAN32", "ZSCAN5B", "ZUFSP", "ZXDA", "ZXDB", "ZXDC", "ZBTB24", "ZNF277", "ZNF512", "MBNL2", "ZC3H8", "ZGPAT", "JRK", "TIGD3", "TIGD4", "TIGD5", "TIGD7", "CAMTA1", "CAMTA2", "CXXC4", "FBXL19", "KDM2A", "TET3", "KMT2B", "DMRTB1", "FLYWCH1", "GATAD2A", "GATAD2B", "ZGLP1", "GRHL3", "GTF2IRD2", "GTF2IRD2B", "HMG20A", "HMGN3", "ADNP", "ADNP2", "LEUTX", "MKX", "NANOGNB", "TPRX1", "ZHX2", "ZHX3", "POU5F2", "HSFX1", "HSFX2", "MSANTD1", "BAZ2B", "MBD3", "MBD4", "MBD6", "PIN1", "SETDB2", "BAZ2A", "MBD1", "MTERF2", "MTERF3", "MTERF4", "DMTF1", "MSANTD4", "MYPOP", "MYSM1", "SNAPC4", "TERB1", "TERF1", "TTF1", "MYRFL", "NFX1", "NFXL1", "AEBP1", "ARHGAP35", "BRF2", "CC2D1A", "CENPA", "CENPS", "CENPT", "CENPX", "CGGBP1", "CHCHD3", "CSRNP1", "CSRNP2", "CSRNP3", "DACH1", "DACH2", "DR1", "DRAP1", "GLMP", "GPBP1", "GPBP1L1", "KCNIP3", "NACC2", "NKRF", "NME2", "PA2G4", "PCGF2", "PCGF6", "PHF19", "PLSCR1", "PREB", "PURB", "PURG", "RAG1", "RBCK1", "REXO4", "SAFB", "SAFB2", "SCMH1", "SKI", "SKIL", "SMYD3", "SNAPC2", "SNAPC5", "SON", "SPEN", "TCF20", "TET2", "THYN1", "TMF1", "TSC22D1", "RFX8", "SP100", "SP110", "SP140", "SP140L", "TBPL1", "THAP10", "THAP11", "THAP2", "THAP3", "THAP4", "THAP5", "THAP6", "THAP7", "THAP8", "THAP9")

         
var <- apply(TF,1,sd)
avg <- apply(TF.count,1,mean)
plot(sort(var^2/avg))
tail(sort(var^2/avg), 20)
variation <- var^2/avg

#to generate a heatmap of pseudotime values
# Define the pseudotime bins
gcdata_CC_regressed_CSS_integrated$pseudotime_bins <- cut(gcdata_CC_regressed_CSS_integrated$slingshot_pseudotime, breaks = 20)
#add bin information to suerat object

#convert to dataframe
tf_list_df <- as.data.frame(tf_list)
pseudotime_bins_df <- as.data.frame(pseudotime_bins)

# Aggregate the gene expression matrix by pseudotime bins
gene_expression_by_bin <- aggregate(tf_list_df, by = list(pseudotime_bins_df), mean)

# Normalize the gene expression values within each gene
gene_expression_norm <- t(scale(t(gene_expression_by_bin)))
# Load the gplots package
library(gplots)

# Plot the heatmap
heatmap.2(gene_expression_norm, col = colorRampPalette(c("blue", "white", "red"))(100), 
          trace = "none", key = FALSE, dendrogram = "none", margins = c(5, 10), 
          main = "Heatmap of scRNA-seq data along pseudotime", xlab = "Pseudotime bins", 
          ylab = "Genes")



#Use the AddModuleScore() function in Seurat to calculate the average expression of your key genes in each bin.
TFs <- read.csv("~/Desktop/all_variable_TFs_AD_scrnaseq.csv", header = TRUE)
# Convert to vector
gene_vector <- as.vector(paste0("'", TFs$Gene, "'", collapse = ","))
print(gene_vector)

#to get gene list####

var_genes <- VariableFeatures(gcdata_CC_regressed_CSS_integrated@assays$RNA)

# Load a data frame of gene symbols and their corresponding transcription factor status
tf_df <- read.csv("~/Desktop/Human_TF.csv", header = TRUE, stringsAsFactors = FALSE)

# Get variable transcription factors
tf_list <- tf_df$Name[tf_df$Name %in% var_genes]

# Print the list of transcription factors
tf_list

# Extract the expression matrix for the selected genes
gene_expression <- as.matrix(gcdata_CC_regressed_CSS_integrated@assays$RNA[ ,tf_list])

gcdata_CC_regressed_CSS_integrated <- AddModuleScore(gcdata_CC_regressed_CSS_integrated, 
                            features = tf_list, 
                            score.name = "gene_expression_TFs", nbins = 21)

#example from chatgpt
gcdata_CC_regressed_CSS_integrated <- AddModuleScore(gcdata_CC_regressed_CSS_integrated, 
                            features = tf_list, 
                            score.name = "gene_expression_TFs", 
                            ctrl = list(bins = gcdata_CC_regressed_CSS_integrated$pseudotime_bins))


DoHeatmap(gcdata_CC_regressed_CSS_integrated, features = tf_list, group.by = "pseudotime_bins")




gcdata_CC_regressed_CSS_integrated <- RunUMAP(gcdata_CC_regressed_CSS_integrated, dims = 1:30)
gcdata_CC_regressed_CSS_integrated <- FindNeighbors(gcdata_CC_regressed_CSS_integrated)
gcdata_CC_regressed_CSS_integrated <- FindClusters(gcdata_CC_regressed_CSS_integrated, resolution = 0.6)
gcdata_CC_regressed_CSS_integrated <- RunUMAP(gcdata_CC_regressed_CSS_integrated, dims = 1:30)
DimPlot(gcdata_CC_regressed_CSS_integrated, reduction = "umap_css", group.by = "seurat_clusters")

#method 2 from chatgpt

# Assuming you have already loaded the Seurat object as "seurat_obj"

# Compute pairwise distances between cells based on pseudotime
pseudotime_dist <- as.dist(1 - cor(t(as.matrix(gcdata_CC_regressed_CSS_integrated@meta.data[, "slingshot_pseudotime"])), use = "pairwise.complete.obs"))


# Cluster cells using fast greedy algorithm
fg <- igraph::fastgreedy.community(as.undirected(pseudotime_dist))

# Divide cells into 20 bins based on cluster assignments
bin_assignments <- paste0("bin", cutree(fg, k = 20))

# Set cell identities in Seurat object based on bin assignments
Idents(gcdata_CC_regressed_CSS_integrated) <- bin_assignments


# Divide dataset into 20 bins based on quantiles of velocity pseudotime
gcdata_CC_regressed_CSS_integrated <- Idents(gcdata_CC_regressed_CSS_integrated) <- paste0("bin", cutree(hclust(dist(gcdata_CC_regressed_CSS_integrated@meta.data$slingshot_pseudotime)), k = 20))

# Compute average gene expression across metacells for each bin
bin_average <- AverageExpression(gcdata_CC_regressed_CSS_integrated, group.by = "ident")

# Compute pairwise Pearson correlation between log-normal gene expression values of each bin
correlation_matrix <- cor(t(log10(bin_average+1)))

# Define distance metric as 1 - correlation coefficient
distance_matrix <- 1 - correlation_matrix

# Perform hierarchical clustering using ward.D2 method
hclust_res <- hclust(as.dist(distance_matrix), method = "ward.D2")

# Manually annotate bins based on resulting clusters
bin_annot <- c(rep("PS cells", 3), "neuroectoderm", "neuroepithelium", "NPCs", "neurons", rep("unknown", 13))
seurat_obj <- SetIdent(seurat_obj, bin_annot[cutree(hclust_res, k = length(unique(bin_annot)))])


VlnPlot(gcdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "timepoint")
VlnPlot(gcdata_mito_filtered_all_10_wk2_15, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "timepoint", pt.size=0 )
#re-filtering genes and mito counts for cells
gcdata_filtered <- subset(gcdata, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .10)

gcdata_mito_filtered_all_10_wk2_15 <- subset(gcdata_mito_filtered, subset = timepoint == "1wk_glia" & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .10 | timepoint == "2wk_glia" & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .15| timepoint == "3wk_glia" & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .10 | timepoint == "4wk_glia" & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .10 | timepoint == "48h_pD" & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .10 | timepoint =="7d_postDox" & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < .10)

gc_1wk_2wk_3wk_4wk <- subset(x = gcdata, subset = (timepoint == "1wk_glia" | timepoint == "2wk_glia" | timepoint == "3wk_glia" | timepoint == "4wk_glia"))
day2_day7 <- subset(x = gcdata, subset = (timepoint == "48h_pD" | timepoint == "7d_postDox"))
dual_smad <- subset(x = gcdata, subset = timepoint == "dualSMAD")
table(gcdata_filtered$timepoint)

wk1_to_wk4_15_mito <- WhichCells(gc_1wk_2wk_3wk_4wk, expression = percent.mito < .15)
day2_day7_10_mito <- WhichCells(day2_day7, expression = percent.mito < .10)
dual_smad_20_mito <- WhichCells(dual_smad, expression = percent.mito < .20)

length(dual_smad_20_mito)

# and subset the object to only keep those cells

gcdata_filtered_wk <- subset(gc_1wk_2wk_3wk_4wk, cells = wk1_to_wk4_15_mito)
gcdata_filtered_day <- subset(day2_day7, cells = day2_day7_10_mito)
gcdata_filtered_dualsmad <- subset(dual_smad, cells = dual_smad_20_mito)

# Merge two Seurat objects
gcdata_mito_filtered <- merge(x = gcdata_filtered_wk, y = gcdata_filtered_day)
gcdata_mito_filtered <- merge(x = gcdata_filtered_wk, y = c(gcdata_filtered_day, gcdata_filtered_dualsmad))

save(gcdata_filtered, file = "~/Desktop/gcdata_filtered_mito_48_day7_10percent_dualSMAD_20_rest_15.RData")



table(gcdata_filtered@meta.data$timepoint)
table(gcdata_filtered_1@meta.data$timepoint)

# plot violins for new data
VlnPlot(gcdata_mito_filtered, features = "percent.mito", group.by = "timepoint")

gcdata_mito_filtered <- subset(gcdata_mito_filtered, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)
DimPlot(gcdata_filtered, reduction = "umap")

#subset seurat object with cells to remove

TFs <- read.csv("~/Desktop/GW20_34_badhuri_cells_to_remove.csv", header = FALSE)

GW20_brain_regions_combined_cells_filtered <- subset(GW20_brain_regions_combined, cells = TFs$V1, invert = TRUE)
save(GW20_brain_regions_combined_cells_filtered, file = "~/Desktop/GW20__brain_regions_combined_cells_filtered_integrated_badhuri.RData")

#differential expression

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
#Positive values indicate that the feature is more highly expressed in the first group.
monocyte.de.markers <- FindMarkers(gcdata_CC_regressed_CSS_integrated, ident.1 = "dualSMAD", ident.2 = "4wk_glia", test.use = "DESeq2")
# view results
head(monocyte.de.markers)
write.csv(monocyte.de.markers, file = "/Users/lindaseong/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/seurat_DE_dualSMAD_4wk_glia_all_cells_deseq2.csv")

gcdata_CC_regressed_CSS_integrated_dualSMAD <- subset(x= gcdata_CC_regressed_CSS_integrated, subset = timepoint == "dualSMAD")

FeaturePlot(object = gcdata_CC_regressed_CSS_integrated, features = c("PAX6", "SOX1", "NES", "SOX9", "FOXG1", "SIX3", "OTX2", "LHX6", "LHX8", "NKX2-1", "SOX6", "SOX2"), cols = c("gray", "red"), pt.size = 0.5, order = TRUE)

FeaturePlot(object = gcdata_CC_regressed_CSS_integrated, features = c("POU5F1", "LIN28A", "DLX1", "DLX2", "GAD1", "GAD2", "SLC32A1", "HES5", "HES6", "SOX2", "MAP2", "SYN1", "MKI67"), cols = c("gray", "red"), pt.size = 0.3, order = TRUE)

g <- DotPlot(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4, features = c("GAL",	"SST", "TAC1",	"TRH",	"NPY",	"ARPP21",	"IGSF1",	"AMIGO2",	"NRGN",	"VIP",	"GRP",	"SCGN",	"GHRH",	"RELN",	"TAC3",	"CARBP1",	"CALY",	"SYT1",	"SYT4",	"CALB1",	"CALB2",	"SEC11C",	"SOX4",	"MEIS2",	"PCSK1N",	"GAP43", "ISL1", "ONECUT1", "ONECUT2", "ONECUT3", "DLX5", "SOX11", "STMN2", "HNRNPH1", "CXADR", "TUBB2A", "DLX1", "RTN1", "STMN1", "TMSB10"), group.by = "seurat_clusters") & theme(axis.text.x=element_text(angle=45, hjust=1, size=9)) & scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu")))

FeaturePlot(object = new_AD7, features = c("GAD1", "GAD2", "SLC32A1", "SLC17A6", "SLC17A7", "SLC6A5", "TH", "CHAT", "SOX2", "TPH2", "NES", "FOXG1", "MAP2", "SYN1", "NKX2-1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)
FeaturePlot(object = new_AD7_mito15, features = c("SIX3", "SYT4", "SLC32A1", "SLC17A6", "MEIS2", "DLX6", "SST", "SOX2", "MAP2", "ARX", "POU2F2", "MEF2C", "LMO4", "TAC1", "ID2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

FeaturePlot(object = gcdata_CC_regressed_CSS_integrated, features = "GRIA3", cols = c("gray", "red"), pt.size = 0.25, order = TRUE) + DimPlot(gcdata_CC_regressed_CSS_integrated, reduction = "umap", group.by = "orig.ident")

t <- subset(new_AD7, subset = FOXG1 > 0)


#attribution plot or stacked bar plot for percent of cells expressing certain genes
map2 <- WhichCells(gcdata_CC_regressed_CSS_integrated, expression = MAP2 > 0)
vgat <- WhichCells(gcdata_CC_regressed_CSS_integrated, expression = SLC32A1 > 0)

gcdata_CC_regressed_CSS_integrated$MAP2_exp<- ifelse(colnames(gcdata_CC_regressed_CSS_integrated) %in% map2, "Pos", "Neg")
gcdata_CC_regressed_CSS_integrated$SLC32A1_exp<- ifelse(colnames(gcdata_CC_regressed_CSS_integrated) %in% vgat, "Pos", "Neg")

gcdata_CC_regressed_CSS_integrated$MAP2_SLC32A1 = with(gcdata_CC_regressed_CSS_integrated,
              ifelse(gcdata_CC_regressed_CSS_integrated$MAP2_exp == "Pos" & gcdata_CC_regressed_CSS_integrated$SLC32A1_exp == "Pos","MAP2+SLC32A1+",
                              ifelse(gcdata_CC_regressed_CSS_integrated$MAP2_exp == "Pos" & gcdata_CC_regressed_CSS_integrated$SLC32A1_exp == "Neg", "MAP2+SLC32A1-",
                                     ifelse(gcdata_CC_regressed_CSS_integrated$MAP2_exp == "Neg" & gcdata_CC_regressed_CSS_integrated$SLC32A1_exp == "Pos", "MAP2-SLC32A1+", "MAP2-SLC32A1-"))))
          
                                         
dittoBarPlot(gcdata_CC_regressed_CSS_integrated, "MAP2_SLC32A1", group.by = "timepoint", scale = "percent", x.reorder = c(4,6,1,2,3,5))

cluster_order <- c("48h_pD", "7d_postDox", "1wk_glia", "2wk_glia", "3wk_glia", "4wk_glia")
gcdata_CC_regressed_CSS_integrated@meta.data$timepoint <- factor(x = gcdata_CC_regressed_CSS_integrated@meta.data$timepoint, levels = cluster_order)

cluster_order <- match(levels(gcdata_CC_regressed_CSS_integrated$timepoint), metaLevels(timepoint, gcdata_CC_regressed_CSS_integrated$timepoint))

then

dittoBarPlot(seurat_object_filtered, var = clusters_to_use, group.by = "your_x_var", var.labels.reorder = var.labels.reorder = order(levels(gcdata_CC_regressed_CSS_integrated@meta.data$timepoint)))

FeaturePlot(object = gcdata_CC_regressed_CSS_integrated, features = c("DDIT3", "EIF2AK3", "XBP1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)
FeaturePlot(object = new_AD7, features = c("DDIT3", "EIF2AK3", "XBP1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

doubles <- subset(gcdata_CC_regressed_CSS_integrated, subset = DDIT3 > 0 & EIF2AK3 > 0)


#regress out ribosomal gene percent
gcdata_CC_regressed_CSS_integrated_stress_0.55_default[["percent.ribo"]] <- PercentageFeatureSet(gcdata_CC_regressed_CSS_integrated_stress_0.55_default, pattern = "^RP[SL]")
#or
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = gcdata_CC_regressed_CSS_integrated), value = TRUE);
percent.rb <- Matrix::colSums(x = GetAssayData(object = gcdata_CC_regressed_CSS_integrated, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = gcdata_CC_regressed_CSS_integrated, slot = 'counts'));
gcdata_CC_regressed_CSS_integrated[['percent.rb']] <- percent.rb;

#remove ribo genes prior to findallmarkers so ribo genes don't show up in DE
genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(seurat_object),
                  value=TRUE, invert=TRUE) #get list of non-ribosomal genes 

# Way2: Doing it manually
ribo_genes <- rownames(gcdata_CC_regressed_CSS_integrated)[grep("^RP[SL]", rownames(gcdata_CC_regressed_CSS_integrated))]
head(ribo_genes, 10)
alldata$percent_ribo <- colSums(alldata@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell

#Then use genes.use in the "features=gene.use"

FeatureScatter(gcdata_CC_regressed_CSS_integrated, "nFeature_RNA", "percent.ribo", group.by = "timepoint", pt.size = 0.15)




library(scvi)
# convert a v3 assay to a v5 assay
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed[["RNA5"]] <- as(object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed[["RNA"]], Class = "Assay5")

#scVI integration
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated <- IntegrateLayers(
  object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, method = scVIIntegration, group.by = "timepoint", orig.reduction = "umap", ndims = 10, nlayers = 2, layers = "counts",
  new.reduction = "integrated.scvi",
  conda_env = "/Users/lindalee/anaconda3/envs/scvi-env", verbose = FALSE
)


gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated <- FindNeighbors(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated, reduction = "integrated.scvi", dims = 1:20)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindClusters(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, resolution = .25)

egcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated <- RunUMAP(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated, reduction = "integrated.scvi", dims = 1:20, reduction.name = "umap.scvi")
DimPlot(
  gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed,
  reduction = "umap",
  group.by = "seurat_clusters")

save(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_scVI_integrated_w_v5.RData")

wrap_plots(c(p1, p2), ncol = 2)


#trying R version of scVI
devtools::install_github("cellgeni/sceasy")
library(reticulate)
library(sceasy)


Sys.setenv(PATH= paste("/Users/lindalee/anaconda3/envs/scvi-env/bin",Sys.getenv()["PATH"],sep=";"))
library(reticulate)
use_condaenv("scvi-env")
py_config()


sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
load("~/Desktop/gcdata_mito_10_wk2_15_filtered_stress_0.56_CC_regressed_20230817.RData")
adata <- convertFormat(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata) # Note generally in Python, dataset conventions are obs x var

# run setup_anndata, use column stim for batch
scvi$model$SCVI$setup_anndata(adata, batch_key = 'timepoint')

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

# get the latent represenation
latent = model$get_latent_representation()

# put it back in our original Seurat object
latent <- as.matrix(latent)
rownames(latent) = colnames(pbmc)
pbmc[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(pbmc))


#load anndata and convert to seurat object
library(Seurat)
devtools::install_github('satijalab/seurat-data')
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratData)
library(SeuratDisk)


d <- GetAssayData(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed[["RNA"]], layer = "counts") 

exp <- FetchData(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, vars = dat$Targets, layer = "data", )

exp <- AverageExpression(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, assays = "RNA", features = dat$Target, group.by = "timepoint", layer = "counts")
head(exp)
write.csv(exp, file = "~/Documents/Ruiqi_scRNAseq/alltimepoints_post_stress_SCENIC_09082023/1wk_postDox_top3_regulonactivity_.7_target_gene_avg_exp.csv")

#seurat V5 integration

mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed <- RunPCA(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed, features = VariableFeatures(object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed), dims = 1:20)
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed, reduction = "pca", group.by = "orig.ident")

#split layers by batch
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed[["RNA"]] <- split(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed[["RNA"]], f = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed$orig.ident)
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed



#harmony
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony <- IntegrateLayers(
  object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony <- RunHarmony(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed, group.by.vars = "orig.ident", reduction = "harmony", assay.use = "RNA", reduction.save = "harmony")

test <- RunUMAP(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, reduction = "harmony", assay = "RNA", dims = 1:20)
test <- FindNeighbors(object = test, reduction = "harmony", dims = 1:20)
test <- FindClusters(test, resolution = .3)
DimPlot(test, reduction = "umap", group.by = "orig.ident")

#fastMNN
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN <- IntegrateLayers(
  object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

#RPCA
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA <- IntegrateLayers(
  object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

#CCA
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_CCA <- IntegrateLayers(
  object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca", group.by = "orig.ident",
  verbose = FALSE
)

mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_SCVI <- IntegrateLayers(
  object = mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, method = scVIIntegration, group.by = "orig.ident", layers = "counts",
  new.reduction = "integrated.scvi",
  conda_env = "/Users/lindalee/anaconda3/envs/scvi-env", verbose = FALSE
)d

DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, reduction = "integrated.mnn")
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, reduction = "integrated.rpca")
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony, reduction = "pca", group.by = "orig.ident")
DimPlot(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, reduction = "pca")


save(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, file = "~/Desktop/AD3_wk4_mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA_integration_20231010.RData")
save(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA, file = "~/Desktop/AD3_wk4_mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA_integration_w_PCA_recalc_prior_to_int_20231005.RData")

mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA <- JoinLayers(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA)
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_harmony_MNN_RPCA


wk4 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "4wk_glia")
wk4_map2_exp <- subset(wk4, subset = MAP2 >0)
wk4_map2_slc32a1_only <- subset(wk4_map2_exp, subset = SLC32A1 >0)
wk4_map2_lc32a1_or_gad <- subset(wk4_map2_exp, subset = SLC32A1 >0 | GAD1 > 0 | GAD2 > 0)



wk4_map2_slc32a1_only <- FindVariableFeatures(wk4_map2_slc32a1_only, nfeatures = 2000)
wk4_map2_slc32a1_only <- RunPCA(wk4_map2_slc32a1_only, features = VariableFeatures(object = wk4_map2_slc32a1_only), dims = 1:20)
#cluster cells
wk4_map2_slc32a1_only <- FindNeighbors(wk4_map2_slc32a1_only, dims = 1:20)
wk4_map2_slc32a1_only <- FindClusters(wk4_map2_slc32a1_only, resolution = .2)
wk4_map2_slc32a1_only <- RunUMAP(wk4_map2_slc32a1_only, dims = 1:20)
DimPlot(wk4_map2_slc32a1_only, reduction = "umap", group.by = "seurat_clusters", pt.size = .3)

save(wk4_map2_slc32a1_only, file = "~/Desktop/wk4_map2_slc32a1_only_mito10_stress_filt_CC_regressed_20231101.RData")


FeaturePlot(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, features = c("POU5F1", "HMG1A", "LIN28A", "PRDX1", "NES", "ID1", "ID3", "SST", "DLX6-AS1", "GNAS"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

FeaturePlot(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, features = c("GNAS", "STMN1", "SYT1", "FOS", "GAP43", "STMN2", "MAP1B", "NEFM", "NEFL"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

FeaturePlot(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, features = c("BASP", "GABARAP", "TUBB3", "HMGA1"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)


gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindClusters(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, resolution = .315)
DimPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed)
DotPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, features = c("POU5F1", "HMGA1", "LIN28A", "PRDX1", "NES", "ID1", "ID3", "SST", "DLX6-AS1", "GNAS", "STMN1", "SYT1", "FOS", "GAP43", "STMN2", "MAP1B", "NEFM", "NEFL", "BASP", "GABARAP", "TUBB3"), group.by = "seurat_clusters") + RotatedAxis() & scale_colour_gradientn(colours =rev(brewer.pal(n = 11, name = "RdBu")))


FeaturePlot(object = wk4_map2_slc32a1_only, features = c("VIP", "SEC11C", "DLX5", "FOXP1", "POU2F2", "ISL1", "MEIS2", "ONECUT2", "ONECUT3", "PBX1", "BCL11B", "FOXA2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

FeaturePlot(object = wk4_map2_slc32a1_only, features = c("SST", "CALB1", "CALB2", "GAD1", "GAD2", "SLC32A1", "SLC17A6", "SLC17A7", "SLC6A5", "TH", "CHAT", "TPH2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)

DotPlot(wk4_interneuron, features = c("SST",	"CALB1", "CALB2",	"AGRP",	"CRH",	"TAC1",	"TRH",	"GAL",	"NPY",	"GHRH",	"GRP",	"VIP",	"SEC11C",	"DLX5",	"FOXP1",	"POU2F2",	"ISL1",	"MEIS2",	"ONECUT2",	"ONECUT3",	"PBX1"), group.by = "seurat_clusters") & theme(axis.text.x=element_text(angle=90, hjust=1, size=9)) & scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu")))

Receptors <- sort(c("VIPR1",	"CRHR2",	"OPRK1",	"SSTR2", "NPY1R",	"NPY5R",	"CRHR1",	"SSTR3",	"OPRL1",	"CNTFR",	"TACR1",	"OXTR",	"OPRD1",	"GPR83",	"GALR1",	"MCHR1",	"NPFFR1",	"NPY2R",	"HCRTR2",	"NMUR1",	"MC4R", "PRLHR", "GHRHR",	"TACR2",	"TACR3",	"GALR2",	"GALR3",	"QRFPR",	"PROKR2",	"VIPR2",	"GRPR",	"KISS1R",	"TRHR",	"PRLR",	"CALCR",	"AR",	"PGR"))

DotPlot(wk4_interneuron, features = Receptors, group.by = "seurat_clusters") & theme(axis.text.x=element_text(angle=90, hjust=1, size=9)) & scale_colour_gradientn(colours =rev(brewer.pal(n = 10, name = "RdBu")))

day2 <- subset(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, subset = timepoint == "48h_pD")

FeaturePlot(object = wk4_map2_slc32a1_only, features = c("NFATC2", "BHLHE41", "MSX1", "ETV1", "IRX5", "LEF1", "TBX2"), cols = c("gray", "red"), pt.size = 0.25, order = TRUE)




