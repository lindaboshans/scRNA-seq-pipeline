# For output from CellRanger < 3.0
data_dir <- '~/Desktop/GW14_Raw_counts/CGE/'
data <- Read10X(data.dir = data_dir)
GW14_CGE = CreateSeuratObject(counts = data, project = "GW14_CGE", min.cells = 3, min.features = 200)
GW14_CGE$region.ident <- "CGE"
GW14_CGE$structure <- "GE"
data_dir <- '~/Desktop/GW14_Raw_counts/LGE/'
data <- Read10X(data.dir = data_dir)
GW14_LGE = CreateSeuratObject(counts = data, project = "GW14_LGE", min.cells = 3, min.features = 200)
GW14_LGE$region.ident <- "LGE"
GW14_LGE$structure <- "GE"
data_dir <- '~/Desktop/GW14_Raw_counts/MGE/'
data <- Read10X(data.dir = data_dir)
GW14_MGE = CreateSeuratObject(counts = data, project = "GW14_MGE", min.cells = 3, min.features = 200)
GW14_MGE$region.ident <- "MGE"
GW14_MGE$structure <- "GE"
data_dir <- '~/Desktop/GW14_Raw_counts/Hypothalamus/'
data <- Read10X(data.dir = data_dir)
GW14_hypothalamus = CreateSeuratObject(counts = data, project = "GW14_hypothalamus", min.cells = 3, min.features = 200)
GW14_hypothalamus$region.ident <- "hypothalamus"
GW14_hypothalamus$structure <- "hypothalamus"
data_dir <- '~/Desktop/GW14_Raw_counts/motor/'
data <- Read10X(data.dir = data_dir)
GW14_motor = CreateSeuratObject(counts = data, project = "GW14_motor", min.cells = 3, min.features = 200)
GW14_motor$region.ident <- "motor"
GW14_motor$structure <- "neocortex"
data_dir <- '~/Desktop/GW14_Raw_counts/Occipital/'
data <- Read10X(data.dir = data_dir)
GW14_V1 = CreateSeuratObject(counts = data, project = "GW14_V1", min.cells = 3, min.features = 200)
GW14_V1$region.ident <- "V1"
GW14_V1$structure <- "neocortex"
data_dir <- '~/Desktop/GW14_Raw_counts/Somato/'
data <- Read10X(data.dir = data_dir)
GW14_somatosensory = CreateSeuratObject(counts = data, project = "GW14_somatosensory", min.cells = 3, min.features = 200)
GW14_somatosensory$region.ident <- "somatosensory"
GW14_somatosensory$structure <- "neocortex"
data_dir <- '~/Desktop/GW14_Raw_counts/Striatum/'
data <- Read10X(data.dir = data_dir)
GW14_striatum = CreateSeuratObject(counts = data, project = "GW14_striatum", min.cells = 3, min.features = 200)
GW14_striatum$region.ident <- "striatum"
GW14_striatum$structure <- "striatum"
data_dir <- '~/Desktop/GW14_Raw_counts/Thalamus/'
data <- Read10X(data.dir = data_dir)
GW14_thalamus = CreateSeuratObject(counts = data, project = "GW14_thalamus", min.cells = 3, min.features = 200)
GW14_thalamus$region.ident <- "thalamus"
GW14_thalamus$structure <- "thalamus"



objects <- ls(pat = "GW14")
object <- objects[-c(1,3,8)]
cell_ids <- c("CGE", "hypo", "LGE", "MGE", "motor", "occipital", "somato", "striatum", "thalamus")

#merge
combined <- merge(GW14_CGE, y = c(GW14_hypothalamus, GW14_LGE, GW14_MGE, GW14_motor, GW14_somatosensory, GW14_striatum, GW14_thalamus, GW14_V1), add.cell.ids = object, project = "GW14_fetal_scRNAseq")
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

save(combined, file = "~/Desktop/GW14_all_regions_merged_20231115.RData")

metadata <- read.csv("~/Desktop/GW14_metadata.csv", header = T, sep = ",")
head(metadata)
dim(metadata)
table(metadata$structure)

#subset based on cluster of specific resolution
#cells.use <- colnames(pbmc_small)[which(pbmc_small[[]]['RNA_snn_res.0.8'] == 1)]
cells.use <- metadata$cell.name
GW14 <- subset(combined, cells = cells.use)
table(GW14$structure)

cell.names <- colnames(GW14)
head(cell.names)
write.csv(cell.names, file = "~/Desktop/GW14_badhuri_cells_to_keep_metadata.csv")

GW14 <- AddMetaData(GW14, metadata = metadata$cell.type, col.name = "celltype")
GW14 <- AddMetaData(GW14, metadata = metadata$clusterv1, col.name = "clusterv1")
GW14 <- AddMetaData(GW14, metadata = metadata$clusterv2...final, col.name = "clusterv2")
head(GW14)

save(GW14, file = "~/Desktop/GW14_all_regions_merged_w_celltype_20231115.RData")

GW14_neurons <- subset(GW14, subset = celltype == "Interneuron" | celltype == "Neuron")
GW14_gad_exp_neurons_and_interneurons <- subset(GW14, subset = GAD1 > 0 & GAD2 > 0)
GW14_gad1_or_2_exp_neurons_and_interneurons <- subset(GW14, subset = GAD1 > 0 | GAD2 > 0)
combined_slc32a1 <- subset(combined, subset = SLC32A1 > 0)
combined_gad1_gad2 <- subset(combined, subset = GAD1 > 0 | GAD2 > 0)
combined_gad1_and_gad2 <- subset(combined, subset = GAD1 > 0 & GAD2 > 0)


#GW20
data_dir <- '~/Desktop/GW20_raw_counts/GW20_hypothalamus/'
data <- Read10X(data.dir = data_dir)
GW20_hypothalamus = CreateSeuratObject(counts = data, project = "GW20_hypothalamus", min.cells = 3, min.features = 200)
GW20_hypothalamus$region.ident <- "hypothalamus"
GW20_hypothalamus$structure <- "hypothalamus"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_caudalthalamus/'
data <- Read10X(data.dir = data_dir)
GW20_caudalthalamus = CreateSeuratObject(counts = data, project = "GW20_caudalthalamus", min.cells = 3, min.features = 200)
GW20_caudalthalamus$region.ident <- "caudalthalamus"
GW20_caudalthalamus$structure <- "thalamus"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_ventral_thalamus/'
data <- Read10X(data.dir = data_dir)
GW20_ventralthalamus = CreateSeuratObject(counts = data, project = "GW20_ventralthalamus", min.cells = 3, min.features = 200)
GW20_ventralthalamus$region.ident <- "ventralthalamus"
GW20_ventralthalamus$structure <- "thalamus"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_dorsalthalamus/'
data <- Read10X(data.dir = data_dir)
GW20_dorsalthalamus = CreateSeuratObject(counts = data, project = "GW20_dorsalthalamus", min.cells = 3, min.features = 200)
GW20_dorsalthalamus$region.ident <- "dorsalthalamus"
GW20_dorsalthalamus$structure <- "thalamus"
data_dir <- '~/Desktop/GW20_raw_counts/GW20V1/'
data <- Read10X(data.dir = data_dir)
GW20_V1 = CreateSeuratObject(counts = data, project = "GW20_V1", min.cells = 3, min.features = 200)
GW20_V1$region.ident <- "V1"
GW20_V1$structure <- "neocortex"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_motor/'
data <- Read10X(data.dir = data_dir)
GW20_motor = CreateSeuratObject(counts = data, project = "GW20_motor", min.cells = 3, min.features = 200)
GW20_motor$region.ident <- "motor"
GW20_motor$structure <- "neocortex"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_somato/'
data <- Read10X(data.dir = data_dir)
GW20_somatosensory = CreateSeuratObject(counts = data, project = "GW20_somatosensory", min.cells = 3, min.features = 200)
GW20_somatosensory$region.ident <- "somatosensory"
GW20_somatosensory$structure <- "neocortex"
data_dir <- '~/Desktop/GW20_raw_counts/GW20PFC/'
data <- Read10X(data.dir = data_dir)
GW20_PFC = CreateSeuratObject(counts = data, project = "GW20_PFC", min.cells = 3, min.features = 200)
GW20_PFC$region.ident <- "PFC"
GW20_PFC$structure <- "neocortex"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_preoptic/'
data <- Read10X(data.dir = data_dir)
GW20_preoptic = CreateSeuratObject(counts = data, project = "GW20_preoptic", min.cells = 3, min.features = 200)
GW20_preoptic$region.ident <- "preoptic"
GW20_preoptic$structure <- "hypothalamus"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_putamen/'
data <- Read10X(data.dir = data_dir)
GW20_putamen = CreateSeuratObject(counts = data, project = "GW20_putamen", min.cells = 3, min.features = 200)
GW20_putamen$region.ident <- "putamen"
GW20_putamen$structure <- "striatum"
data_dir <- '~/Desktop/GW20_raw_counts/GW20CGE/'
data <- Read10X(data.dir = data_dir)
GW20_CGE = CreateSeuratObject(counts = data, project = "GW20_CGE", min.cells = 3, min.features = 200)
GW20_CGE$region.ident <- "CGE"
GW20_CGE$structure <- "GE"
data_dir <- '~/Desktop/GW20_raw_counts/GW20LGE/'
data <- Read10X(data.dir = data_dir)
GW20_LGE = CreateSeuratObject(counts = data, project = "GW20_LGE", min.cells = 3, min.features = 200)
GW20_LGE$region.ident <- "LGE"
GW20_LGE$structure <- "GE"
data_dir <- '~/Desktop/GW20_raw_counts/GW20MGE/'
data <- Read10X(data.dir = data_dir)
GW20_MGE = CreateSeuratObject(counts = data, project = "GW20_MGE", min.cells = 3, min.features = 200)
GW20_MGE$region.ident <- "MGE"
GW20_MGE$structure <- "GE"
data_dir <- '~/Desktop/GW20_raw_counts/GW20_claustrum/'
data <- Read10X(data.dir = data_dir)
GW20_claustrum = CreateSeuratObject(counts = data, project = "GW20_claustrum", min.cells = 3, min.features = 200)
GW20_claustrum$region.ident <- "claustrum"
GW20_claustrum$structure <- "claustrum"

objects <- ls(pat = "GW20_")
object <- objects[-c(1:8)]
cell_ids <- c("caudate", "CGE", "dorsalthalamus", "hypo", "LGE", "MGE", "motor", "NucAcc", "parietal", "parVZ", "PFC", "PFCVZ", "putamen", "V1", "ventralthalamus")

#merge
GW20_combined <- merge(GW20_caudalthalamus, y = c(GW20_CGE, GW20_claustrum, GW20_dorsalthalamus, GW20_hypothalamus, GW20_LGE, GW20_MGE, GW20_motor, GW20_PFC, GW20_preoptic, GW20_putamen, GW20_somatosensory, GW20_V1, GW20_ventralthalamus), add.cell.ids = objects, project = "GW20")
GW20_combined <- NormalizeData(GW20_combined)
GW20_combined <- FindVariableFeatures(GW20_combined, selection.method = "vst", nfeatures = 2000)

GW20_combined <- ScaleData(GW20_combined)
GW20_combined <- RunPCA(GW20_combined, npcs = 30)
GW20_combined <- RunUMAP(GW20_combined, reduction = "pca", dims = 1:30)
GW20_combined <- FindNeighbors(GW20_combined, reduction = "pca", dims = 1:30)
GW20_combined <- FindClusters(GW20_combined, resolution = 0.5)

DimPlot(GW20_combined)
save(GW20_combined, file = "~/Desktop/GW20_new_individual_merged_regions_20231117.RData")
load("~/Desktop/GW20_new_individual_merged_regions_20231114.RData")

GW20_meta <- read.csv("~/Desktop/GW20_metadata.csv", header = TRUE, sep = ",")
dim(GW20_meta)
head(GW20_meta)


#subset based on cluster of specific resolution
#cells.use <- colnames(pbmc_small)[which(pbmc_small[[]]['RNA_snn_res.0.8'] == 1)]
cells.use <- GW20_meta$cell.name
GW20 <- subset(GW20_combined, cells = cells.use)
table(GW20$structure)
meta <- GW20@meta.data
dim(meta)
write.csv(meta, file = "~/Desktop/test.csv")
meta <- read.csv(file = "~/Desktop/test.csv", header = TRUE, sep = ",")
head(meta)

GW20 <- AddMetaData(GW20, metadata = GW20_meta$cell.type, col.name = "celltype")
GW20 <- AddMetaData(GW20, metadata = GW20_meta$clusterv1, col.name = "clusterv1")
GW20 <- AddMetaData(GW20, metadata = GW20_meta$clusterv2...final, col.name = "clusterv2")
head(GW20@meta.data)

save(GW20, file = "~/Desktop/GW20_new_individual_merged_regions_w_celltype_20231114.RData")


GW20_neurons <- subset(GW20, subset = celltype == "Interneuron" | celltype == "Neuron")
GW20_gad1_or_2_exp_neurons_and_interneurons <- subset(GW20_neurons, subset = GAD1 > 0 | GAD2 > 0)
GW20_gad_exp_neurons_and_interneurons <- subset(GW20_neurons, subset = GAD1 > 0 & GAD2 > 0)
GW20_slc32a1 <- subset(GW20_neurons, subset = SLC32A1 > 0)
GW20_gad1_gad2 <- subset(GW20_neurons, subset = GAD1 > 0 | GAD2 > 0)
GW20_gad1_and_gad2 <- subset(GW20_neurons, subset = GAD1 > 0 & GAD2 > 0)
GW20_gad_or_slc32a1 <- subset(GW20_neurons, subset = GAD1 > 0 | GAD2 > 0 | SLC32A1 > 0)
#subset for only areas in GW20_34 dataset
GW20_gad_exp_neurons_and_interneurons_areas_removed <- subset(GW20_gad_exp_neurons_and_interneurons, subset = region.ident == "claustrum" | region.ident == "somatosensory", invert = TRUE)

save(GW20, file = "~/Desktop/GW20_new_individual_merged_regions_w_celltype_neurons_only_20231117.RData")




#GW20_31
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_hypo/'
data <- Read10X(data.dir = data_dir)
GW20_31_hypothalamus = CreateSeuratObject(counts = data, project = "GW20_31_hypothalamus", min.cells = 3, min.features = 200)
GW20_31_hypothalamus$region.ident <- "hypothalamus"
GW20_31_hypothalamus$structure <- "hypothalamus"
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_V1/'
data <- Read10X(data.dir = data_dir)
GW20_31_V1 = CreateSeuratObject(counts = data, project = "GW20_31_V1", min.cells = 3, min.features = 200)
GW20_31_V1$region.ident <- "V1"
GW20_31_V1$structure <- "neocortex"
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_parietal/'
data <- Read10X(data.dir = data_dir)
GW20_31_parietal = CreateSeuratObject(counts = data, project = "GW20_31_parietal", min.cells = 3, min.features = 200)
GW20_31_parietal$region.ident <- "parietal"
GW20_31_parietal$structure <- "neocortex"
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_CGE/'
data <- Read10X(data.dir = data_dir)
GW20_31_CGE = CreateSeuratObject(counts = data, project = "GW20_31_CGE", min.cells = 3, min.features = 200)
GW20_31_CGE$region.ident <- "CGE"
GW20_31_CGE$structure <- "GE"
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_LGE/'
data <- Read10X(data.dir = data_dir)
GW20_31_LGE = CreateSeuratObject(counts = data, project = "GW20_31_LGE", min.cells = 3, min.features = 200)
GW20_31_LGE$region.ident <- "LGE"
GW20_31_LGE$structure <- "GE"
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_MGE/'
data <- Read10X(data.dir = data_dir)
GW20_31_MGE = CreateSeuratObject(counts = data, project = "GW20_31_MGE", min.cells = 3, min.features = 200)
GW20_31_MGE$region.ident <- "MGE"
GW20_31_MGE$structure <- "GE"
data_dir <- '~/Desktop/GW20_31_raw_counts/GW20_31_PFC/'
data <- Read10X(data.dir = data_dir)
GW20_31_PFC = CreateSeuratObject(counts = data, project = "GW20_31_PFC", min.cells = 3, min.features = 200)
GW20_31_PFC$region.ident <- "PFC"
GW20_31_PFC$structure <- "neocortex"


objects <- ls(pat = "GW20_31")
object <- objects[-c(1:8)]
cell_ids <- c("caudate", "CGE", "dorsalthalamus", "hypo", "LGE", "MGE", "motor", "NucAcc", "parietal", "parVZ", "PFC", "PFCVZ", "putamen", "V1", "ventralthalamus")

#merge
GW20_31_combined <- merge(GW20_31_CGE, y = c(GW20_31_hypothalamus, GW20_31_LGE, GW20_31_MGE, GW20_31_parietal, GW20_31_PFC, GW20_31_V1), add.cell.ids = objects, project = "GW20_31")

GW20_31_combined <- NormalizeData(GW20_31_combined)
GW20_31_combined <- FindVariableFeatures(GW20_31_combined, selection.method = "vst", nfeatures = 2000)
GW20_31_combined <- ScaleData(GW20_31_combined)
GW20_31_combined <- RunPCA(GW20_31_combined, npcs = 30)
GW20_31_combined <- RunUMAP(GW20_31_combined, reduction = "pca", dims = 1:30)
GW20_31_combined <- FindNeighbors(GW20_31_combined, reduction = "pca", dims = 1:30)
GW20_31_combined <- FindClusters(GW20_31_combined, resolution = 0.5)

DimPlot(GW20_31_combined)
save(GW20_31_combined, file = "~/Desktop/GW20_31_individual_merged_regions_w_celltype_20231117.RData")
load("~/Desktop/GW20_new_individual_merged_regions_20231114.RData")

GW20_meta <- read.csv("~/Desktop/GW20_31_metadata.csv", header = TRUE, sep = ",")
dim(GW20_meta)
head(GW20_meta)
table(GW20_meta$structure)


#subset based on cluster of specific resolution
#cells.use <- colnames(pbmc_small)[which(pbmc_small[[]]['RNA_snn_res.0.8'] == 1)]
cells.use <- GW20_meta$cell.name
GW20_31 <- subset(GW20_31_combined, cells = cells.use)
table(GW20_31$structure)
meta <- GW20_31@meta.data
dim(meta)
write.csv(meta, file = "~/Desktop/test.csv")
meta <- read.csv(file = "~/Desktop/test.csv", header = TRUE, sep = ",")
head(meta)

GW20_31 <- AddMetaData(GW20_31, metadata = meta$cell.type, col.name = "celltype")
GW20_31 <- AddMetaData(GW20_31, metadata = meta$clusterv1, col.name = "clusterv1")
GW20_31 <- AddMetaData(GW20_31, metadata = meta$clusterv2...final, col.name = "clusterv2")
head(GW20_31@meta.data)

save(GW20_31, file = "~/Desktop/GW20_new_individual_merged_regions_w_celltype_20231114.RData")


GW20_31_neurons <- subset(GW20_31, subset = celltype == "Interneuron" | celltype == "Neuron")
GW20_gad1_or_2_exp_neurons_and_interneurons <- subset(GW20_31_neurons, subset = GAD1 > 0 | GAD2 > 0)
GW20_31_neurons_gad_exp_neurons_and_interneurons <- subset(GW20_31_neurons, subset = GAD1 > 0 & GAD2 > 0)
GW20_slc32a1 <- subset(GW20_neurons, subset = SLC32A1 > 0)
GW20_gad1_gad2 <- subset(GW20_neurons, subset = GAD1 > 0 | GAD2 > 0)
GW20_gad1_and_gad2 <- subset(GW20_neurons, subset = GAD1 > 0 & GAD2 > 0)
GW20_gad_or_slc32a1 <- subset(GW20_neurons, subset = GAD1 > 0 | GAD2 > 0 | SLC32A1 > 0)
#subset for only areas in GW20_34 dataset
GW20_gad_exp_neurons_and_interneurons_areas_removed <- subset(GW20_gad_exp_neurons_and_interneurons, subset = region.ident == "claustrum" | region.ident == "somatosensory", invert = TRUE)

save(GW20_31_neurons, file = "~/Desktop/GW20_31_merged_regions_w_celltype_20231117_neurons_only_20231117.RData")

#chnage headings to include age for PCA plto

