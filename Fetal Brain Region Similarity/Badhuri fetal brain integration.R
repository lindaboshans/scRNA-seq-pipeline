#this script creates seurat objects of each region, merges by timepoint, and integrates as needed. 


#function for automating this process

library(Seurat)
library(SeuratDisk)
library(dplyr)


automateSeuratObjectCreation <- function(timepoint, base_dir) {
  raw_counts_dir <- file.path(base_dir, paste0(timepoint, '_raw_counts'))
  
  structure_mapping <- c(
    "MGE" = "GE", 
    "CGE" = "GE", 
    "LGE" = "GE",
    "hypothalamus" = "hypothalamus",
    "thalamus" = "thalamus",
    "somato" = "neocortex",
    "temporal" = "neocortex",
    "parietal" = "neocortex",
    "occipital" = "neocortex",
    "motor" = "neocortex",
    "V1" = "neocortex",
    "wholePFC" = "neocortex",
    "PFC_CP" = "neocortex",
    "V1_CP" = "neocortex",
    "motor_CP" = "neocortex",
    "PFC" = "neocortex",
    "striatum" = "striatum",
    "nucleusaccumbens" = "striatum",
    "caudate" = "striatum",
    "putamen" = "striatum",
    "ventral_thalamus" = "thalamus",
    "dorsalthalamus" = "thalamus"
  )
  
  region_folders <- list.dirs(raw_counts_dir, full.names = TRUE, recursive = FALSE)
  
  seurat_objects <- list()
  for (folder in region_folders) {
    region_name <- basename(folder)
    specific_region <- gsub(paste0(timepoint, '_'), '', region_name)
    data <- Read10X(data.dir = folder)
    
    seurat_object <- CreateSeuratObject(counts = data, project = paste0(timepoint, '_', specific_region), min.cells = 3, min.features = 200)
    seurat_object$region.ident <- specific_region
    seurat_object$structure <- unname(structure_mapping[specific_region])
    
    seurat_objects[[specific_region]] <- seurat_object
  }
  
  # Merge all Seurat objects into one
  combined_object <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = names(seurat_objects), project = paste0(timepoint, "_combined"))
  
  combined_object[["RNA"]] <- JoinLayers(combined_object[["RNA"]])
  
  combined_object <- NormalizeData(combined_object) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData() %>% 
    RunPCA(npcs = 30) %>% 
    RunUMAP(reduction = "pca", dims = 1:30) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(resolution = 0.5)
  
  saveRDS(combined_object, file = file.path(base_dir, paste0(timepoint, '_new_individual_merged_regions_', Sys.Date(), '.RData')))
  
  return(combined_object)
}


# usage:
GW18_combined <- automateSeuratObjectCreation(timepoint = "GW18", base_dir = "~/Desktop")



#subset for cell ids provided by Badhuri paper with metadata info 
addMetadataToSeurat <- function(seurat_object, metadata_file_path, timepoint) {
  # Read the metadata CSV file
  metadata <- read.csv(metadata_file_path, header = TRUE, sep = ",")
  
  # Filter metadata for the specific timepoint if required
  # Assuming there's a column in your metadata that specifies the timepoint
  metadata_subset <- metadata[metadata$individual == timepoint, ]
  
  # Adjust cell names in metadata to match the format of cell names in the Seurat object
  #metadata_subset$cell.name <- gsub("^GW\\d+_\\d+_", "", metadata_subset$cell.name) #<- this code is for "GWXX_X" or "GWXX_XX" formatted datasets
  metadata_subset$cell.name <- gsub("^GW\\d{2}(_\\d+)?_", "", metadata_subset$cell.name) #<- this code is for "GWXX" formatted datasets 
  metadata_subset$cell.name <- paste0(metadata_subset$cell.name, "-1")
  
  # Subset Seurat object to keep only the cells present in the metadata
  cells_use <- intersect(colnames(seurat_object), metadata_subset$cell.name)
  seurat_object <- subset(seurat_object, cells = cells_use)
  
  # Ensure the metadata is in the same order as the cells in the Seurat object
  metadata_subset <- metadata_subset[match(cells_use, metadata_subset$cell.name), ]
  
  # Add metadata columns to the Seurat object
  seurat_object <- AddMetaData(seurat_object, metadata = metadata_subset$cell.type, col.name = "celltype")
  seurat_object <- AddMetaData(seurat_object, metadata = metadata_subset$clusterv1, col.name = "clusterv1")
  seurat_object <- AddMetaData(seurat_object, metadata = metadata_subset$`clusterv2...final`, col.name = "clusterv2")
  
  # Return the updated Seurat object
  return(seurat_object)
}

# Example usage

GW18_combined <- addMetadataToSeurat(GW18_combined, "~/Desktop/metadata_badhuri.csv", "GW18")

# Save the updated Seurat object
saveRDS(GW18_combined, file = "~/Desktop/GW18_all_regions_merged_w_celltype_20240305.rds")




#subset for gad exp neurons
GW18_combined_neurons <- subset(GW18_combined, subset = celltype == "Interneuron" | celltype == "Neuron")
GW18_combined_gad_exp_neurons_and_interneurons <- subset(GW18_combined_neurons, subset = GAD1 > 0 | GAD2 > 0)
saveRDS(GW18_combined_gad_exp_neurons_and_interneurons, file = "~/Desktop/GW18_all_regions_merged_w_celltype_neurons_interneurons_gadexp_20240305.rds")


load("~/Desktop/GW20_31_merged_regions_w_celltype_20231117_neurons_only_20231117.RData")
GW20_31_combined_gad_exp_neurons_and_interneurons <- subset(GW20_31_neurons, subset = GAD1 > 0 | GAD2 > 0)

load("~/Desktop/GW20_new_individual_merged_regions_w_celltype_neurons_only_20231117.RData")
GW20_gad_exp_neurons_and_interneurons <- subset(GW20_neurons, subset = GAD1 > 0 | GAD2 > 0)

GW18_2_combined_gad_exp_neurons_and_interneurons <- readRDS("~/Desktop/GW18_2_all_regions_merged_w_celltype_neurons_interneurons_gadexp_20240305.rds")
GW18_combined_gad_exp_neurons_and_interneurons <- readRDS("~/Desktop/GW18_all_regions_merged_w_celltype_neurons_interneurons_gadexp_20240305.rds")
load("~/Desktop/GW20_31_merged_regions_w_celltype_20231117_neurons_only_20231117.RData")
GW20_31_neurons_gad_exp_neurons_and_interneurons <- subset(GW20_31_neurons, subset = GAD1 > 0 | GAD2 > 0)
GW19_neurons_gad_exp_neurons_and_interneurons<- readRDS("~/Desktop/GW19_all_regions_merged_w_celltype_neurons_interneurons_gadexp_20240305.rds")
load("~/Documents/Ruiqi_scRNAseq/RData_files/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_neurons_only_badhuri_060802023.RData")
GW20_brain_regions_combined_cells_filtered_neurons_gad_exp <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = GAD1 > 0 | GAD2 > 0)


#merge all timepoints into one seurat object
cell.ids <- c("GW18", "GW18_2", "GW18", "GW20", "GW20_31", "GW20_34")
all_fetal_ages_gad_exp_neurons_interneurons <- merge(x = GW18_combined_gad_exp_neurons_and_interneurons, y = c( GW18_2_combined_gad_exp_neurons_and_interneurons, GW19_neurons_gad_exp_neurons_and_interneurons, GW20_gad_exp_neurons_and_interneurons, GW20_31_neurons_gad_exp_neurons_and_interneurons, GW20_brain_regions_combined_cells_filtered_neurons_gad_exp), add.cell.ids = cell.ids, project = "Badhuri_fetal_regions")

all_fetal_ages_gad_exp_neurons_interneurons[["RNA"]] <- JoinLayers(all_fetal_ages_gad_exp_neurons_interneurons[["RNA"]])

all_fetal_ages_gad_exp_neurons_interneurons <- subset(all_fetal_ages_gad_exp_neurons_interneurons, subset = structure == "claustrum", invert = TRUE)

saveRDS(all_fetal_ages_gad_exp_neurons_interneurons, "~/Desktop/badhuri_GW18_19_20_fetal_brain_regions_merged_gad_exp_neurons_interneurons_only_20240306.rds")








#integration of fetal brain data with AD iN data
#############

#integrate brain region datasets and object of interest for PCA analysis
#chnage headings of seurat objects prior to merge to include age for PCA plot

# Access and modify the column
metadata <- GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons@meta.data
structure <- metadata$structure
structure_collapsed <- metadata$structure_collapsed
test <- paste0("GW20_34_", structure)
collapsed <- paste0("GW20_", structure)
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- AddMetaData(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, metadata = test, col.name = "structure")
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- AddMetaData(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, metadata = collapsed, col.name = "structure_collapsed")
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- AddMetaData(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, metadata = "GW20_34", col.name = "dataset")
head(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons)

metadata <- GW20_new_indvidual_neurons_GADexp_neurons_interneurons@meta.data
structure <- metadata$structure
test <- paste0("GW20_", structure)
collapsed <- paste0("GW20_", structure)
GW20_new_indvidual_neurons_GADexp_neurons_interneurons <- AddMetaData(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, metadata = test, col.name = "structure")
GW20_new_indvidual_neurons_GADexp_neurons_interneurons <- AddMetaData(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, metadata = collapsed, col.name = "structure_collapsed")
GW20_new_indvidual_neurons_GADexp_neurons_interneurons <- AddMetaData(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, metadata = "GW20", col.name = "dataset")
head(GW20_new_indvidual_neurons_GADexp_neurons_interneurons)

metadata <- AD7_map2_slc32a1@meta.data
structure <- metadata$seurat_clusters
test <- paste0("PGP1_AD7", structure)
collapsed <- paste0("PGP1_AD7", structure)
AD7_map2_slc32a1 <- AddMetaData(AD7_map2_slc32a1, metadata = structure, col.name = "structure")
AD7_map2_slc32a1 <- AddMetaData(AD7_map2_slc32a1, metadata = collapsed, col.name = "structure_collapsed")
AD7_map2_slc32a1 <- AddMetaData(AD7_map2_slc32a1, metadata = "PGP1_AD7", col.name = "dataset")
head(AD7_map2_slc32a1)

#merge
combined <- merge(AD7_map2_slc32a1, y = c(GW20_new_indvidual_neurons_GADexp_neurons_interneurons, GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons), add.cell.ids = c("PGP1_AD7", "GW20_34", "GW20"),  project = "brain_similarity")

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

saveRDS(combined, file = "~/Desktop/PGP1_AD7_merged_w_GW20_GW20_34_fetal_brain_regions_20240229.rds")



combined <- readRDS("~/Desktop/PGP1_AD7_merged_w_GW20_GW20_34_fetal_brain_regions_20240229.rds")
#integration 

# split the merged  dataset into a list  seurat objects 
ifnb.list <- SplitObject(combined, split.by = "dataset")


# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list, dims = 1:30)

#perform integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
PGP1_AD7_fetal_brain_seurat_integrated <- IntegrateData(anchorset = immune.anchors)

saveRDS(PGP1_AD7_fetal_brain_seurat_integrated, "~/Desktop/PGP1_AD7_GW20_fetal_brain_seurat_integrated_immune_anchors_20240229.rds")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(PGP1_AD7_fetal_brain_seurat_integrated) <- "integrated"

# re-join layers after integration
PGP1_AD7_fetal_brain_seurat_integrated[["RNA"]] <- JoinLayers(PGP1_AD7_fetal_brain_seurat_integrated[["RNA"]])

# Run the standard workflow for visualization and clustering
Seurat.list <- "PGP1_AD7_fetal_brain_seurat_integrated"

for (i in 1:length(Seurat.list)) {
  Seurat.object.name <- Seurat.list[i]
  Seurat.object <- get(Seurat.object.name)
  Seurat.object <- ScaleData(Seurat.object, verbose = FALSE)
  Seurat.object <- RunPCA(Seurat.object, npcs = 30, dims = 1:30)
  Seurat.object <- RunUMAP(Seurat.object, reduction = "pca", dims = 1:30)
  Seurat.object <- RunUMAP(Seurat.object, dims = 1:30)
  Seurat.object <- FindNeighbors(Seurat.object)
  Seurat.object <- FindClusters(Seurat.object, resolution = 0.4)
}

for (i in 1:length(Seurat.list)) {
  Seurat.object.name <- Seurat.list[i]
  Seurat.object <- get(Seurat.object.name)
  Seurat.object <- ScaleData(Seurat.object, verbose = FALSE)
  Seurat.object <- RunPCA(Seurat.object, dims = 1:20)
  Seurat.object <- RunUMAP(Seurat.object, reduction = "pca", dims = 1:20)
  Seurat.object <- RunUMAP(Seurat.object, dims = 1:20)
  Seurat.object <- FindNeighbors(Seurat.object)
  Seurat.object <- FindClusters(Seurat.object, resolution = 0.4)
}

saveRDS(Seurat.object, "~/Desktop/PGP1_AD7_GW20_fetal_brain_seurat_integrated_immune_anchors_20240229.rds")

PGP1_AD7_fetal_brain_seurat_integrated <- readRDS("~/Desktop/PGP1_AD7_GW20_fetal_brain_seurat_integrated_immune_anchors_20240229.rds")

# Visualization
DimPlot(PGP1_AD7_fetal_brain_seurat_integrated, reduction = "pca", group.by = "structure")
DimPlot(PGP1_AD7_fetal_brain_seurat_integrated, reduction = "umap", group.by = "structure", pt.size = .001)
p2 <- DimPlot(gcdata_seurat_integrated, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


DimPlot(object = PGP1_AD7_fetal_brain_seurat_integrated, dims = c(2, 5), reduction = "pca", group.by = "structure")



