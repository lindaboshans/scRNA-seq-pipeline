if (!require(scTypeR)) {
  install_github("grisslab/scTypeR")}

devtools::install_url("http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
                      args="--configure-args='--with-jags-include=/Users/casallas/homebrew/opt/jags/include/JAGS
--with-jags-lib=/Users/casallas/homebrew/opt/jags/lib'
"
)


BiocManager::install("devtools")
library(devtools)

BiocManager::install("infercnv")
BiocManager::install("ComplexHeatmap")
install.packages("fastqcr")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install the scTyper package
devtools::install_github("omicsCore/scTyper")
install.packages("processx")

library(scTypeR)
library("scTyper") 

# phenotype data for test data
pheno.fn = system.file("~/Desktop/pheno_info_public.csv", package = "scTyper")
pheno.fn = read.csv("~/Desktop/sctyper.csv", stringsAsFactors = FALSE)
pheno.fn = "~/Desktop/sctyper.csv"
read.csv(~/Desk)

# user defined marker list
marker=list(T_cell=c("CD2", "CD3D", "CD3E", "CD3G"),
            Fibroblast=c("FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1"),
            Macrophage=c("CD14", "CD163", "CD68", "FCGR2A", "CSF1R"),
            Dendritic_cell=c("CD40", "CD80", "CD83", "CCR7"),
            Mast_cell=c("CMA1", "MS4A2", "TPSAB1", "TPSB2"),
            Myocyte=c("ACTA1", "ACTN2", "MYL2", "MYH2"),
            Endothelial.cell=c("PECAM1", "VWF", "ENG"),
            B_Plasma_cell=c("SLAMF7", "CD79A", "BLNK", "FCRL5"),
            Epithelial_cell=c("KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8",
                              "KRT16", "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3","EPCAM", "SFN"))

# single study
marker="Puram.2017.HNSCC.TME"
# multiple study
marker=c("Kawai.2018.Liver", "Li.2017.CRC.TME")

#Users can also set cell markers by combining identifiers of interest. Following is examples:
# Identifier
  marker=c("Costea.2013.OSCC.CAF:Normal cell:NF", "Costea.2013.OSCC.CAF:Cancer cell:CAF_D",
           "Elyada.2019.PDAC.CAF:Cancer cell:iCAF", "Elyada.2019.PDAC.CAF:Cancer cell:myCAF")

marker=c("Molnar.2003.Dorsolateral prefrontal cortex:Normal cell:Glutamatergic neuron", 
         "Leong.2016.Brain:Normal cell:Neuron", 
         "Company.abcam.Undefined:Normal cell:GABAergic neuron", "Yousef.2018.Brain:Normal cell:Neuron", "Moroni.2007.Brain:Normal cell:Neuron", 
         "Company.abcam.Undefined:Normal cell:Immature neuron", "Company.abcam.Undefined:Normal cell:Serotonergic neuron", "Tagliafierro.2016.Brain:Normal cell:Neuron",
         "Chu.2017.Undefined:Normal cell:Neuronal progenitor cell", "Zhong.2018.Embryonic prefrontal cortex:Normal cell:Interneuron", 
        "Okazaki.2019.Inferior colliculus:Normal cell:Neural stem cell", "Haus.2015.Brain:Normal cell:Neural stem cell", 
         "Sandberg.2016.Brain:Normal cell:Neural stem cell", "Sandberg.2016.Brain:Normal cell:Neural progenitor cell")



marker=c("Company.abcam.Undefined:Normal cell:GABAergic neuron", "Zhong.2018.Embryonic prefrontal cortex:Normal cell:Interneuron")


D1 <- read.csv("~/Desktop/D2 MSN genes link stri paper.csv", header = FALSE)
D2 <- read.csv("~/Desktop/mature MSN DE genes.csv")
D3 <- read.csv("~/Desktop/MGE interneuron genes linc paper .csv")
D4 <- read.csv("~/Desktop/Migrating interneuron DE genes linc paper.csv")
D5 <- read.csv("~/Desktop/ventral CGE int DE genes linc paper.csv")
hyp <- read.csv("~/Desktop/hypo neuron genes biorxiv paper .csv")

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
                                                    "INA", "LDHA", "LHX6", "MAF", "MEF2C", "NEUROD1", "NEUROD6", "NRXN3", "PDE4DIP", "PPP1R17", "PTPRZ1", "VIM"),
                                      D2_MSN=D1,
                                      D1_MSN_mature=D2,
                                      MGE_int=D3,
                                      Migrating_int=D4,
                                      vCGE_int=D5,
                                      Hypothalamus=hyp)

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
                                                    "INA", "LDHA", "LHX6", "MAF", "MEF2C", "NEUROD1", "NEUROD6", "NRXN3", "PDE4DIP", "PPP1R17", "PTPRZ1", "VIM"),
                                      Striatal_interneuron=genes$Striatum,
                                      Hypothalamus=hyp)


#run sctyper
celltyped.seurat = scTyper(seurat.object=gcdata, marker=marker,
                         output.name = "scTyper_day7PD.output",
                         pheno.fn = pheno.fn,
                         qc = FALSE,
                         run.cellranger=FALSE,
                         norm.seurat=FALSE,
                         cell.typing.method="Average",
                         level="cluster",
                         run.inferCNV=FALSE,
                         proj.name = "scTyper", report.mode=TRUE, mc.cores = 2)

day7pd <- cell.typing.seurat(gcdata, marker=marker_cellmarker_and_pangloDB, cell.typing.method = "NTP", level = "cluster", wd, slot = "data", assay = "RNA",
                             ntp.dir = "~/Desktop",
                             NTP.gene.filter.cutoff = 0.3,
                             NTP.distance = c("cosine", "correlation"),
                             NTP.norm.method = "none",
                             mc.cores = 1) 


expr = GetAssayData(gcdata, slot = "data")
markerList=get.markerList(marker_cellmarker_and_pangloDB)
markerList=lapply(markerList, function(a) intersect(a, rownames(gcdata)))
marker.mean=sapply(markerList, function(a) as.matrix(if(length(unlist(a))==1){expr[a,]}else{colMeans(expr[a,], na.rm=T)}))
cell.type=apply(marker.mean, 1, function(a) names(sort(-a)[1]))
gcdata$average=cell.type
df=table(cell.type, gcdata$seurat_clusters)
rownames(df)
c.name=apply(df,2, function(a) rownames(df)[order(-a)[1]])
c.cell.type=sapply(gcdata$seurat_clusters, function(a) c.name[a])
gcdata$cell.type=as.factor(c.cell.type)

DimPlot(gcdata, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cell.type')

expr = GetAssayData(day2, slot = "data")
markerList=get.markerList(marker_cellmarker_and_pangloDB)
markerList=lapply(markerList, function(a) intersect(a, rownames(day2)))
marker.mean=sapply(markerList, function(a) as.matrix(if(length(unlist(a))==1){expr[a,]}else{colMeans(expr[a,], na.rm=T)}))
cell.type=apply(marker.mean, 1, function(a) names(sort(-a)[1]))
day2$average=cell.type
df=table(cell.type, day2$seurat_clusters)
rownames(df)
c.name=apply(df,2, function(a) rownames(df)[order(-a)[1]])
c.cell.type=sapply(day2$seurat_clusters, function(a) c.name[a])
day2$cell.type.2.sctypeR=as.factor(c.cell.type)

DimPlot(day2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cell.type.2.sctypeR')
