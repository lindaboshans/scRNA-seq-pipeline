# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = c("Embryo", "Brain") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

marker_cellmarker_and_pangloDB = list(ESC=c("LEFTY2", "MYC", "FOXD3", "CD24", "CDH1", "PRDM5", "ANO6", "CDK8", "SLC46A2", "TEX19", "ZSCAN10", "KCNIP3", "PCGF2", "GAL", "PITX2", "BCL3", "SALL4", "DPPA5A", "TDGF1", "HOXB5", "SUMO2", "BNIP3", "GJC1", "TRIM6", "SOX15", "LMNA", "FBXO15", "EPCAM", "ITGB1", "CD59", "PECAM1", "PROM1", "KLF4", "ZFP42", "TAF8", "STAT3", "HMGA2", "NACC1", "NR6A1", "LEF1", "L1TD1", "KHDC3L", "GDF3", "DPPA4", "SMAD1", "SMAD5", "SMAD9", "SMAD4", "SMAD2", "SMAD3", "CTNNB1", "HES1", "DNMT3B", "TRIM28", "NANOG", "POU5F1", "ZFP42", "SOX2", "TERT", "PODXL", "FUT4", "LIN28A", "CHST10", "NGFR", "EPCAM", "CDH3", "CD9"),
                                      GABAergic_neuron=c("SLC6A1", "GABBR2", "GAD1", "GADD45B", "PAX2", "SLC32A1", "GAD2", "VIP", "SST", "PVALB", "TNFAIP8L3", "SEMA3C", "MYBPC1", "PARM1", "DLX1", "DLX2", "DLX5", "DLX6"),
                                      Cholinergic_neuron=c("CHAT", "SLC18A3", "ACHE", "SLC5A7", "TAC1", "BRCA1", "ACLY"),
                                      Neuroblast=c("NTNG1", "TSHZ1", "SP9", "SP8", "PROK2", "ASCL1", "DLL3", "SALL3", "DRAXIN", "ERBB4", "TRIM32", "IGFBPL1", "EOMES", "PBX1", "CUX2", "ISL1", "PHF1", "NEUROG1", "NEUROG2", "NKX6-1", "MAPK7", "PROKR2", "POU6F2", "ITGA6", "EFNA2", "DCX", "GRM5", "EZH2", "DLX2", "MARK2", "LRP8", "VLDLR", "DAB1", "NCAN", "EPHA4", "TACC3", "EGFR", "TUBB3", "SIRT1", "NCAM1", "NEUROD1", "SCGN", "KLK3", "NES", "WT1", "ZBED4", "PROS1", "CNR1"),
                                      Immature_neuron=c("PAX6", "PCLO", "CUX2", "SLC12A5", "POMC", "STIM2", "RAX", "CUL5", "HIVEP2", "FUT9", "TBR1", "PROX1", "NTRK2", "RND2", "FEZF2", "SOX5", "ROBO1", "EPHA7", "EPHA3", "YY2", "SLC1A1", "NOTCH3", "NOTCH1", "LHX6", "NES", "RBFOX3", "CREB1", "NEUROD1", "ASCL1", "NTRK1", "SYN1", "NF1", "BCL2", "SLC6A5", "FOXD2", "ICAM5", "SSTR2", "NEUROG2"),
                                      Neural_stem_cell=c("NGFR", "SOX1", "SOX2", "SOX9", "PROM1", "NGFR", "SLITRK6", "NES", "ADGRG1", "DCX", "NEUROD1", "NCAM1", "DPYSL3", "VIM", "NOTCH1", "MSI1", "MKI67", "CXCR4", "ASCL1"),
                                      Neural_progenitor_cell=c("NGFR", "GFAP", "NES", "SOX1", "SOX2", "SOX9", "PROM1", "NGFR", "SLITRK6", "NES", "ADGRG1", "DCX", "NEUROD1", "NCAM1", "DPYSL3", "ATAD2", "ANP32E", "APOLD1", "AURKB", "BTG3", "BUB1", "BUB3", "CCNA2", "CCNB1", "CCNB2", "CDC20", "CDCA3", "CDK1", "CDKN3", "CKAP2", "EMX1", "EMX2", "EOMES", "FBXO5", "FEN1", "GAP43", "GLI3", "HES1", "HSPB1", "ID4", "KIF22", "LMO3", "LMO4", "MAPT", "MCM2", "MCM3", "MAD2L1", "MEF2C", "MELK", "MFAP2", "MYT1L", "NRXN1", 
                                                               "NOTCH1", "PAX6", "OTX1", "PTTG1", "RCN1", "SYT1", "SOX3", "TGIF1", "TPX2", "TUBA1B", "TUBB4B", "ZEB2", "NEUROG1", "NEUROG2", "DLX2", "DLX5", "SP9", "POU3F4", "LHX9", "SHOX2", "DBX1", "ARX", "HOPX", "HES5", "ASCL2", "S100A6"),
                                      Interneuron=c("ARX", "CELF4", "CKS2", "CXCR4", "DLX2", "DLX1",
                                                    "DLX5", "ERBB4", "FDFT1", "GAD1", "GAD2", "GNG5","HMGB2", "IGFBPL1",
                                                    "INA", "LDHA", "LHX6", "MAF", "MEF2C", "NEUROD1", "NEUROD6", "NRXN3", "PDE4DIP", "PPP1R17", "PTPRZ1", "VIM"))

marker_cellmarker_and_pangloDB = list(GABAergic_neuron=c("SLC6A1", "GABBR2", "GAD1", "GADD45B", "PAX2", "SLC32A1", "GAD2", "VIP", "SST", "PVALB", "TNFAIP8L3", "SEMA3C", "MYBPC1", "PARM1", "DLX1", "DLX2", "DLX5", "DLX6"),
                                      Immature_neuron=c("PAX6", "PCLO", "CUX2", "SLC12A5", "POMC", "STIM2", "RAX", "CUL5", "HIVEP2", "FUT9", "TBR1", "PROX1", "NTRK2", "RND2", "FEZF2", "SOX5", "ROBO1", "EPHA7", "EPHA3", "YY2", "SLC1A1", "NOTCH3", "NOTCH1", "LHX6", "NES", "RBFOX3", "CREB1", "NEUROD1", "ASCL1", "NTRK1", "SYN1", "NF1", "BCL2", "SLC6A5", "FOXD2", "ICAM5", "SSTR2", "NEUROG2"),
                                      Interneuron=c("ARX", "CELF4", "CKS2", "CXCR4", "DLX2", "DLX1",
                                                    "DLX5", "ERBB4", "FDFT1", "GAD1", "GAD2", "GNG5","HMGB2", "IGFBPL1",
                                                    "INA", "LDHA", "LHX6", "MAF", "MEF2C", "NEUROD1", "NEUROD6", "NRXN3", "PDE4DIP", "PPP1R17", "PTPRZ1", "VIM"))

marker_cellmarker_and_pangloDB = list(SST_neuron=c("GAD1", "GAD2", "LHX6", "SST"), #CR+?
                                      PV_neuron=c("GAD1", "GAD2", "LHX6", "PVALB"),
                                      VIP_neuron=c("GAD1", "GAD2", "ADRAB2", "VIP"),
                                      non_VIP_neuron=c("GAD1", "GAD2", "ADRAB2", "LAMP5")) #incl reeln? 


# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = gcdata[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

es.max = sctype_score(scRNAseqData = gcdata[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = marker_cellmarker_and_pangloDB) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(gcdata@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(gcdata@meta.data[gcdata@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(gcdata@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/2] = "Unknown"
print(sctype_scores[,1:3])

gcdata@meta.data$scType_w_custom_markers_10cluster = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  gcdata@meta.data$scType_w_custom_markers10[gcdata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(gcdata, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'scType_w_custom_markers10')    



r#to get bubbles graph

# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
nodes <- nodes[-71, ]

mygraph <- graph_from_data_frame(edges)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

ncells = sum(gcdata@meta.data$seurat_clusters==cl)

scater::multiplot(DimPlot(gcdata, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)
