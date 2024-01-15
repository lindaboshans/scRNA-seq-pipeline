install.packages("devtools")
devtools::install_github("pcahan1/singleCellNet")
library(singleCellNet)

DimPlot(gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2, reduction = "umap", group.by = "timepoint")
load("/Users/lindaseong/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_neurons_only_badhuri_060802023.RData")
load("/Users/lindalee/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/gcdata_with_dualSMAD_mito_10_and_15_and20_CC_regressed_CSS_integrated_1000_union_genes_res_0.9_w_slingshot_pseudotime_20230601.RData")
metadata <- read.csv("~/Documents/Ruiqi_scRNAseq/GW20_34_badhuri_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
df <- as.data.frame(GW20_brain_regions_combined_cells_filtered@meta.data)
write.csv(df, "~/Desktop/GW20_34_metadata_all_cells.csv", row.names = TRUE)
load("~/Documents/Ruiqi_scRNAseq/RData_files/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_neurons_only_badhuri_060802023.RData")
load("/Users/lindalee/Library/Mobile Documents/com~apple~CloudDocs/Desktop/gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4_w_annotation.RData")

metadata_nodup <- metadata[!duplicated(metadata$cells),]

#subset seurat object with cells to remove

TFs <- read.csv("/Users/lindaseong/Desktop/cells_to_remove.csv", header = FALSE)

GW20_brain_regions_combined_cells_filtered <- subset(GW20_brain_regions_combined, cells = TFs$V1, invert = TRUE)
save(GW20_brain_regions_combined_cells_filtered, file = "~/Desktop/GW20_brain_regions_combined_cells_filtered_integrated_badhuri_060802023.RData")

GW20_brain_regions_combined_cells_filtered <- AddMetaData(object = GW20_brain_regions_combined_cells_filtered, metadata = metadata_nodup$clusterv2...final, col.name = 'clusterv2')
save(GW20_brain_regions_combined_cells_filtered, file = "~/Desktop/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_badhuri_060802023.RData")

GW20_brain_regions_combined_cells_filtered_neurons <- subset(x = GW20_brain_regions_combined_cells_filtered, subset = (celltype == "Neuron" | celltype == "Interneuron" ))
save(GW20_brain_regions_combined_cells_filtered_neurons, file = "~/Desktop/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_neurons_only_badhuri_060802023.RData")
load("~/Documents/Ruiqi_scRNAseq/RData_files/GW20_brain_regions_combined_cells_filtered_integrated_w_celltype_neurons_only_badhuri_060802023.RData")

GW20_brain_regions_combined_cells_filtered_neurons_GADexp_int_nonGADexp_neurons_small_regions_removed <- subset(x = GW20_brain_regions_combined_cells_filtered_neurons_GADexp_int_nonGADexp_neurons, subset = (region.ident == "hypothalamus" | region.ident == "nucAcc" | region.ident == "putamen" | region.ident == "ventralthalamus"), invert = TRUE)

GW20_brain_regions_combined_cells_filtered_interneurons <- subset(x= GW20_brain_regions_combined_cells_filtered_neurons, subset = celltype == "Interneuron")
GW20_brain_regions_combined_cells_filtered_neurons <- subset(x= GW20_brain_regions_combined_cells_filtered_neurons, subset = celltype == "Neuron")
GW20_brain_regions_combined_cells_filtered_interneurons_slc32a1_only <- subset (x = GW20_brain_regions_combined_cells_filtered_neurons, subset = SLC32A1 > 0)
GW20_brain_regions_combined_cells_filtered_neurons_GADexp_int_nonGADexp_neurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = celltype == "Interneuron" & GAD1 > 0 & GAD2 > 0 | celltype == "Neuron" & GAD1 == 0 & GAD2 == 0)
GW20_brain_regions_combined_cells_filtered_neurons_nonGADexp_neurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = celltype == "Neuron" & GAD1 >= 0 & GAD2 >= 0, invert = TRUE)
GW20_brain_regions_combined_cells_filtered_neurons_GADexp_intneurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = celltype == "Interneuron" & GAD1 <= 0 & GAD2 <= 0, invert = TRUE)
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_all_interneurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = celltype == "Neuron" & GAD1 <= 0 & GAD2 <= 0, invert = TRUE)
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = GAD1 > 0 & GAD2 > 0)
GW20_brain_regions_combined_cells_filtered_GAD1_or_2_exp_neurons_interneurons <- subset(GW20_brain_regions_combined_cells_filtered_neurons, subset = GAD1 > 0 | GAD2 > 0)
GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons_regions_removed <- subset(GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons, subset = region.ident == "PFCVZ" | region.ident == "parVZ", invert = TRUE)
GW20_new_indiv_GADexp_neurons_all_interneurons <- subset(GW20_gad1_and_gad2, subset = celltype == "Neuron" & GAD1 <= 0 & GAD2 <= 0, invert = TRUE)
GW20_new_and_34_integrated <- subset(wk4_fetal_brain_seurat_integrated, subset = dataset == "GW20" | dataset == "GW20_34")
t <- as.data.frame(wk4_fetal_brain_seurat_integrated$structure)
write.csv(t, file = "~/Desktop/test.csv")
t <- read.csv(file = "~/Desktop/test.csv", header = TRUE, sep = ",")
wk4_fetal_brain_seurat_integrated <- AddMetaData(wk4_fetal_brain_seurat_integrated, metadata = t$wk4, col.name = "structure_collapsed")
save(wk4_fetal_brain_seurat_integrated, file = "~/Desktop/wk4_fetal_brain_seurat_integrated_gad_exp_only_20231121.RData")


GW20_brain_regions_combined <- NormalizeData(GW20_brain_regions_combined)
GW20_brain_regions_combined <- FindVariableFeatures(GW20_brain_regions_combined, nfeatures = 2000)

GW20_brain_regions_combined <- ScaleDatas(GW20_brain_regions_combined, verbose = FALSE)
GW20_brain_regions_combined <- RunPCA(GW20_brain_regions_combined, npcs = 30, verbose = FALSE)
GW20_brain_regions_combined <- RunUMAP(GW20_brain_regions_combined, dims = 1:30)
GW20_brain_regions_combined <- FindNeighbors(GW20_brain_regions_combined, dims = 1:30)
GW20_brain_regions_combined <- FindClusters(GW20_brain_regions_combined, resolution = 0.5)

cells <- as.data.frame(GW20_brain_regions_combined_cells_filtered_neurons_GADexp_int_nonGADexp_neurons$orig.ident)
write.csv(cells, file='~/Desktop/GW20_brain_regions_integrated_badhuri_cellnames.csv', quote=F, row.names=F)

#integrate seurat object to SCN analysis
#exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
#training data 
GW20_new_and_34_integrated$cell <- colnames(GW20_new_and_34_integrated)

GW20_brain_regions_combined_counts = extractSeurat(GW20_new_and_34_integrated, exp_slot_name = "counts")
stTM = GW20_brain_regions_combined_counts$sampTab
expTMraw = GW20_new_and_34_integrated@assays$RNA@counts

#querydata
#conservative interneuron parameters
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only <- subset(gcdata_CC_regressed_CSS_integrated, subset = SLC32A1 > 0 )
gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2 <- subset(gcdata_CC_regressed_CSS_integrated, subset = SLC32A1 > 0 & GAD1 > 0 & GAD2 > 0)
gcdata_CC_regressed_CSS_integrated_slc32a1_or_gad1_or_gad2 <- subset(gcdata_CC_regressed_CSS_integrated, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)
t_filtered_iN_genes <- subset(t, subset = SLC32A1 > 0 & GAD1 > 0 & GAD2 > 0)
DimPlot(new_AD7_slc32a1_or_gad1_or_gad2, reduction = "umap", group.by = "seurat_clusters")

new_AD7_slc32a1_or_gad1_or_gad2 <- subset(new_AD7, subset = SLC32A1 > 0 | GAD1 > 0 | GAD2 > 0)
new_AD7_slc32a1_and_gad1_and_gad2 <- subset(new_AD7, subset = SLC32A1 > 0 & GAD1 > 0 & GAD2 > 0)

wk4_map2_slc32a1_only$cell <- colnames(wk4_map2_slc32a1_only)
gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2_counts = extractSeurat(wk4_map2_slc32a1_only, exp_slot_name = "counts")
stPark = gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2_counts$sampTab
expPark = wk4_map2_slc32a1_only@assays$RNA@counts

t <- readRDS(file = "~/Downloads/GSE233574_OrganoidScreen_processed_SeuratObject")


#example data

#download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda", "sampTab_Park_MouseKidney_062118.rda")
#download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expMatrix_Park_MouseKidney_Oct_12_2018.rda", "expMatrix_Park_MouseKidney_Oct_12_2018.rda")
#download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expMatrix_TM_Raw_Oct_12_2018.rda", "expMatrix_TM_Raw_Oct_12_2018.rda")
#download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda", "sampTab_TM_053018.rda")

#stPark = utils_loadObject("sampTab_Park_MouseKidney_062118.rda")
#expPark = utils_loadObject("expMatrix_Park_MouseKidney_Oct_12_2018.rda")
#expTMraw = utils_loadObject("~/Downloads/stList_Darminis.rda")
#stTM = utils_loadObject("sampTab_TM_053018.rda")
genesPark = rownames(expPark)
stTM<-droplevels(stTM)
commonGenes = intersect(rownames(expTMraw), genesPark)
length(commonGenes)
expTMraw = expTMraw[commonGenes,]
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=stTM, ncells=100, dLevel="structure_collapsed")
stTrain = stList[[1]]
#commonCells = intersect(colnames(expTMraw), rownames(stTrain))
#expTrain = expTMraw[,commonCells]
expTrain = expTMraw[,rownames(stTrain)]

#train the classifier
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 25, nRand = 70, nTrees = 1000, nTopGenePairs = 70, dLevel = "structure_collapsed", colName_samp = "cell"))

#apply to held out data
#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="structure_collapsed") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = expTMraw[commonGenes,rownames(stTest)]

#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 70)

#assess classifier
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "structure_collapsed", classQuery = "structure_collapsed", nRand = 70)

plot_PRs(tm_heldoutassessment)
plot_metrics(tm_heldoutassessment)
tm_heldoutassessment$AUPRC_w
tm_heldoutassessment$AUPRC_wc
tm_heldoutassessment$accuracy
#these are the training results for dataset GW20_brain_regions_combined_cells_filtered_GADexp_neurons_interneurons
#for dataset GW20_brain_regions_combined_cells_filtered_neurons(or unfiltered nueron/interneuron), ncells = 700 is optimal
#ncells = 50 gives accuracy of 0.909375, w of 0.9399497 and wc of 0.9436289
#ncells = 100 gives accuracy of 0.9245614, w of 0.9624134 and wc of 0.9605896
#ncells = 200 gives accuracy of 0.8560794, w of 0.9170037 and wc of 0.9230584
#ncells = 400 gives accuracy of 0.884252, w of 0.9657005 and wc of 0.9569897
#ncells = 500 gives accuracy of 0.8707006, w of 0.9597485 and wc of 0.9490099
#ncells = 600 gives accuracy of 0.8807487, w of 0.966006 and wc of 0.9569068
#ncells = 700 gives accuracy of 0.8879566, w of 0.9656978 and wc of 0.9591429
#ncells = 800 gives accuracy of 0.901223, w of 0.9590662 and wc of 0.9586846
#ncells = 900 gives accuracy of 0.9043127, w of 0.9539226 and wc of 0.9639678
#ncells = 1000 gives accuracy of 0.9079966, w of 0.9426253 and wc of 0.9687906
#ncells = 1200 gives accuracy of 0.922407, w of 0.8741558 and wc of 0.9828869


#classification result heatmap

#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names

nrand = 50
sla = as.vector(stTest$structure)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created

sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)

#attribution plot
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="structure", sid="cell")

#visualize average top pairs gene expression for training data
gpTab = compareGenePairs(query_exp = expTest, training_exp = expTrain, training_st = stTrain, classCol = "structure", sampleCol = "cell", RF_classifier = class_info$cnProc$classifier, numPairs = 20, trainingOnly= TRUE)

train = findAvgLabel(gpTab = gpTab, stTrain = stTrain, dLevel = "structure")

hm_gpa_sel(gpTab, genes = class_info$cnProc$xpairs, grps = train, maxPerGrp = 50)



##QUERY##
#apply to Park et al query data
#expPark = utils_loadObject("expMatrix_Park_MouseKidney_Oct_12_2018.rda") 

nqRand = 70
system.time(crParkall<-scn_predict(class_info[['cnProc']], expPark, nrand=nqRand))

plot_attr(crParkall, stPark, nrand=nqRand, sid="cell", dLevel="seurat_clusters")

visualization
sgrp = as.vector(stPark$structure)
names(sgrp) = as.vector(stPark$cell)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

# heatmap classification result
sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)

#classification annotation assignment
# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.

stPark <- get_cate(classRes = crParkall, sampTab = stPark, dLevel = "seurat_clusters", sid = "cell", nrand = nqRand)

#classificaion result violin plot
sc_violinClass(sampTab = stPark, classRes = crParkall, sid = "cell", dLevel = "seurat_clusters", addRand = nqRand)
sc_violinClass(stPark, crParkall, sid = "cell", dLevel = "timepoint", ncol = 12, sub_cluster = "4wk_glia")

#skyline plot of classification results
library(viridis)
stKid2 = addRandToSampTab(crParkall, stPark, "timepoint", "cell")
skylineClass(crParkall, "dualSMAD", stKid2, "timepoint",.25, "cell")


system.time(umPrep_HS<-prep_umap_class(crParkall, stPark, nrand=nqRand, dLevel="seurat_clusters", sid="cell", topPC=5))
plot_umap(umPrep_HS)

gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2_week4 <- subset(gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2, subset = timepoint == "4wk_glia")
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4 <- RunPCA(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4)
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4 <- RunUMAP(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4, reduction = "pca", dims = 1:20)
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4 <- RunUMAP(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4, assay = "RNA", dims = 1:20)
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4 <- FindNeighbors(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4, dims = 1:20)
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4 <- FindClusters(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4, dims = 1:20, resolution = 0.25)
DimPlot(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4, reduction = "umap", group.by = "seurat_clusters")
g <- FindAllMarkers(gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2_week4)
g

gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only <- subset(gcdata_CC_regressed_CSS_integrated, subset = SLC32A1 > 0 )
gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only_week4 <- subset(gcdata_CC_regressed_CSS_integrated_slc32a1_exp_only, subset = timepoint == "4wk_glia")

DimPlot(gcdata_CC_regressed_CSS_integrated, reduction = "umap", group.by = "seurat_clusters")

organoid_interneurons <- subset(organoid_neurons, subset = class2 == "neuron.gaba")
f <- sum(as.numeric(WhichCells(object = organoid_interneurons, expression = ONECUT3 > 0)))

FeaturePlot(stri_11pcw, features = c("ONECUT1", "ONECUT2", "ONECUT3"), reduction = "umap") + DimPlot(stri_11pcw, reduction = "umap", group.by = "Cell_type")

DimPlot(stri_11pcw, reduction = "umap", group.by = "Cell_type") 