library(Seurat)

AD7_map2_slc32a1 <- readRDS("~/Desktop/PGP1_WT2_AD7_map2_and_slc32a1_only_mito5_CC_regressed_SCT_20250228.rds")
all_fetal_ages_gad_exp_neurons_interneurons <- readRDS("~/Desktop/badhuri_GW18_19_20_fetal_brain_regions_merged_gad_exp_neurons_interneurons_only_SCT_20240306.rds")
load("~/Desktop/wk4_map2_slc32a1_only_mito10_stress_filt_CC_regressed_w_pyscenicAUC_GW20_annotation_20231114.RData")
organoids <- readRDS("~/Desktop/GSE233574_OrganoidScreen_processed_SeuratObject.rds")
AD_combined_map2_slc32a1 <- readRDS("~/Desktop/PGP1_AD_methods_map2_slc32a1_MNN_int_CC_ribo_hb_regressed_SCT_w_ref_map_fetal_mid_hind_hypo_20240319.rds")
load("~/Desktop/RData files/AD3_mito20_stress_filtered_CC_and_ribo_regressed_20230928.RData")

organoids <- subset(organoids, subset = region == "non-specific", invert = TRUE)
organoids_neurons <- subset(organoids, subset = class == "neuron")
organoid_interneurons <- subset(organoids_neurons, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)

AD3_mito20_stress_filt_CC_and_ribo_regressed <- subset(AD3_mito20_stress_filt_CC_and_ribo_regressed, subset = MAP2 > 0 & SLC32A1 > 0)
AD3_mito20_stress_filt_CC_and_ribo_regressed <- SCTransform(AD3_mito20_stress_filt_CC_and_ribo_regressed, verbose = FALSE)
organoid_interneurons <- SCTransform(organoid_interneurons, verbose = FALSE)
all_fetal_ages_gad_exp_neurons_interneurons_nohypo <- subset(all_fetal_ages_gad_exp_neurons_interneurons, subset = structure == "hypothalamus", invert =TRUE)
table(all_fetal_ages_gad_exp_neurons_interneurons_nothalamus$structure)

AD3_mito20_stress_filt_CC_and_ribo_regressed <- RunPCA(AD3_mito20_stress_filt_CC_and_ribo_regressed, verbose = FALSE)
AD3_mito20_stress_filt_CC_and_ribo_regressed <- RunUMAP(AD3_mito20_stress_filt_CC_and_ribo_regressed, dims = 1:20, verbose = FALSE)
AD3_mito20_stress_filt_CC_and_ribo_regressed <- FindNeighbors(AD3_mito20_stress_filt_CC_and_ribo_regressed, dims = 1:20, verbose = FALSE)
AD3_mito20_stress_filt_CC_and_ribo_regressed <- FindClusters(AD3_mito20_stress_filt_CC_and_ribo_regressed, verbose = FALSE, resolution = 0.1)
DimPlot(AD_combined_map2_slc32a1, label = TRUE, reduction = "umap")

organoid_interneurons <- RunPCA(organoid_interneurons, verbose = FALSE)
all_fetal_ages_gad_exp_neurons_interneurons <- RunUMAP(all_fetal_ages_gad_exp_neurons_interneurons, dims = 1:30, verbose = FALSE, return.model=TRUE)
organoid_interneurons <- FindNeighbors(organoid_interneurons)
organoid_interneurons <- FindClusters(organoid_interneurons, verbose = FALSE)
DimPlot(all_fetal_ages_gad_exp_neurons_interneurons, label = TRUE)

saveRDS(AD7_map2_slc32a1, "~/Desktop/PGP1_WT2_AD7_map2_and_slc32a1_only_mito5_CC_regressed_SCT_20250228.rds")
saveRDS(all_fetal_ages_gad_exp_neurons_interneurons, "~/Desktop/badhuri_GW18_19_20_fetal_brain_regions_merged_gad_exp_neurons_interneurons_only_SCT_20240306.rds")
saveRDS(wk4_map2_slc32a1_only, "~/Desktop/wk4_map2_slc32a1_only_mito10_stress_filt_CC_regressed_w_pyscenicAUC_GW20_annotation_SCT_20240306.rds")

DefaultAssay(all_fetal_ages_gad_exp_neurons_interneurons_nothalamus) <- "SCT"
DefaultAssay(AD7_map2_slc32a1) <- "SCT"
DefaultAssay(wk4_map2_slc32a1_only) <- "SCT"


#reference mapping using SCT
anchors <- FindTransferAnchors(
  reference = all_fetal_ages_gad_exp_neurons_interneurons,
  query = AD3_mito20_stress_filt_CC_and_ribo_regressed,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)

AD3_mito20_stress_filt_CC_and_ribo_regressed <- MapQuery(
  anchorset = anchors,
  query = AD3_mito20_stress_filt_CC_and_ribo_regressed,
  reference = all_fetal_ages_gad_exp_neurons_interneurons,
  refdata = list(
    structure = "structure",
    age = "orig.ident"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)


AD methods map2 slc32a1 reference mapping w fetal herb hypo midbrain on umap

DimPlot(AD_combined_map2_slc32a1, group.by = "predicted.structure", reduction = "umap.mnn") + DimPlot(AD_combined_map2_slc32a1, reduction = "umap.mnn", group.by = "RNA_snn_res.0.23")
DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, group.by = "predicted.structure", reduction = "umap") + DimPlot(AD3_mito20_stress_filt_CC_and_ribo_regressed, reduction = "umap", group.by = "RNA_snn_res.0.2")
DimPlot(AD_combined_map2_slc32a1, group.by = "predicted.age", reduction = "umap.mnn")

table(wk4_map2_slc32a1_only$predicted.structure)
table(AD3_mito20_stress_filt_CC_and_ribo_regressed$predicted.structure)
saveRDS(AD_combined_map2_slc32a1, "~/Desktop/PGP1_AD_methods_map2_slc32a1_MNN_int_CC_ribo_hb_regressed_SCT_w_ref_map_fetal_mid_hind_hypo_20240320.rds")
saveRDS(AD7_map2_slc32a1, "~/Desktop/PGP1_WT2_AD7_map2_and_slc32a1_only_mito5_CC_regressed_SCT_w_refmap_fetal_mid_hind_hypo_20240320.rds")
saveRDS(AD3_mito20_stress_filt_CC_and_ribo_regressed, "~/Desktop/AD3_mito20_stress_filtered_CC_and_ribo_regressed_SCT_w_refmap_fetal_mid_hind_hypo_20240320.rds")

DefaultAssay(AD_methods) <- "RNA"
AD_methods <- FindNeighbors(AD_methods, dims = 1:20, verbose = FALSE)
AD_methods <- FindClusters(AD_methods, verbose = FALSE, resolution = 0.3)
DimPlot(AD_methods, reduction = "umap.mnn")

DefaultAssay(all_fetal_ages_gad_exp_neurons_interneurons) <- "RNA"
DefaultAssay(AD7_map2_slc32a1) <- "RNA"

#reference mapping using lognormalize
anchors <- FindTransferAnchors(
  reference = all_fetal_ages_gad_exp_neurons_interneurons,
  query = wk4_map2_slc32a1_only,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:20
)

wk4_map2_slc32a1_only <- MapQuery(
  anchorset = anchors,
  query = wk4_map2_slc32a1_only,
  reference = all_fetal_ages_gad_exp_neurons_interneurons,
  refdata = list(
    celltype.log = "celltype",
    structure.log = "structure",
    age.log = "orig.ident"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(wk4_map2_slc32a1_only, group.by = "predicted.structure.log", reduction = "umap")
DimPlot(wk4_map2_slc32a1_only, group.by = "predicted.structure", reduction = "umap")



#stacked bar plot for results
# Extract metadata
metadata <- AD_combined_map2_slc32a1@meta.data

# Count the number of cells for each condition and cell type
data_to_plot <- metadata %>%
  group_by(predicted.structure, orig.ident) %>%
  summarise(cell_count = n()) %>%
  ungroup() %>%
  mutate(proportion = cell_count / sum(cell_count))

# Calculate the proportions
data_to_plot <- data_to_plot %>%
  group_by(predicted.structure) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup()

ggplot(data_to_plot, aes(x = predicted.structure, y = proportion, fill = orig.ident)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Brain Region", y = "Proportion", fill = "AD Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
