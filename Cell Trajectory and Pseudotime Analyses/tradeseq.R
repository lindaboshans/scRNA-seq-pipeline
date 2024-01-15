if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
BiocManager::install("tradeSeq")

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
load("~/Desktop/gcdata_percent_mito_10_and_15_CC_regressed_CSS_integrated_union_top1000genes_res_0.9_w_slingshot_pseudotime_20230518.RData")

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))

counts_matrix <- GetAssayData(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, assay='RNA', slot='counts')
counts <- as.matrix(counts_matrix)
df_counts <- as.data.frame(counts)
rm(counts_matrix)
sce <- as.SingleCellExperiment(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed,  assay = "RNA")
sce <- slingshot(sce, clusterLabels = 'timepoint', reducedDim = 'PHATE')
crv <- SlingshotDataSet(sce)
data(celltype, package = "tradeSeq")

#fit negative binomial
set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:8, 
                   nGenes = 200, verbose = T)
#fitGAM
set.seed(7)
BiocParallel::register(BiocParallel::SerialParam())
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 2 # use 2 cores
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- FindVariableFeatures(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, nfeatures = 10000)
genes <- as.vector(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@assays$RNA@var.features)
pseudotime <- slingPseudotime(crv, na = FALSE)
pseudotime <- as.matrix(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$phate_slingshot_pseudotime1)
cellWeights <- slingCurveWeights(crv)
conditions <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$timepoint
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights, 
              nknots = 8, verbose = TRUE, genes = genes, parallel=TRUE)

#below is more efficient
sce <- fitGAM(counts = counts,
              sds = crv, nknots = 8, genes = genes, verbose = TRUE)

table(rowData(sce)$tradeSeq$converged)

saveRDS(sce, "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_fitGAM_sce_object_10kvargenes_09122023.rds")

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce)$assocRes
write.csv(assocRes, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_10kvargenes_assocRes_09122023.csv")

saveRDS(sce, "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_fitGAM_sce_object_10kvargenes_w_assocRes.rds")
sce <- readRDS("~/Desktop/gcdata_mito_10and15_CC_regressed_CSS_integrated_0.9_fitGAM_sce_object_10kvargenes_w_assocRes_startRes_05252023.rds")

library(tradeSeq)
assocRes <- rowData(sce)$assocRes
mockGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05)
]


length(mockGenes)
write.csv(mockGenes, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_10kvargenes_assocRes_sig_genes_only_09122023.csv")

#visualization of DE genes
### based on mean smoother
library(pheatmap)
yhatSmooth <- predictSmooth(sce, gene = mockGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)

yhatSmooth <- predictSmooth(sce, gene = mockGenes, nPoints = 50)

## the hierarchical trees constructed here, can also be used for 
## clustering of the genes according to their average expression pattern.
cl <- sort(cutree(heatSmooth$tree_row, k = 6))
table(cl)

conditions <- colData(sce)$pheno$treatment_id
pt1 <- colData(sce1)$slingshot$pseudotime

### based on fitted values (plotting takes a while to run)
yhatCell <- predictCells(sce, gene=mockGenes)
df <- as.data.frame(yhatCell)
yhatCellMock <- yhatCell[,conditions == "Mock"]

# order according to pseudotime
ooMock <- order(pseudotime, decreasing=FALSE)
ooMock <- order(pseudotime[conditions == "Mock"], decreasing=FALSE)
yhatCell <- yhatCell[,ooMock]
pheatmap(t(scale(t(yhatCell))), cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames=FALSE)

#discovering progenitor marker genes
startRes <- startVsEndTest(sce)
rowData(sce)$startRes <- startRes
sig_genes_startvsend <-  rownames(startRes)[
  which(p.adjust(startRes$pvalue, "fdr") <= 0.05)
]

sig_genes_startvsend_all_info <- 
length(sig_genes_startvsend)

#startRes <- rowData(sce)$startRes
#sig_genes_p_adjust <- p.adjust(startRes$pvalue, method = "fdr")
#write.csv(startRes, file = "~/Desktop/tradeseq_mito_10_10kvargenes_startRes_sig_genes_only.csv")
#write.csv(sig_genes_p_adjust, file = "~/Desktop/tradeseq_mito_10_10kvargenes_startRes_sig_genes_only_adj_values_only.csv")

write.csv(startRes, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_10kvargenes_startvsEndRes_05252023.csv")
write.csv(sig_genes_startvsend, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_10kvargenes_startvsEndRes_sig_genes_only.csv")

saveRDS(sce, "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_fitGAM_sce_object_10kvargenes_w_pseudotime_scaled_0_to_1_10262023.rds")
sce <- readRDS("~/Documents/Ruiqi_scRNAseq/tradeseq_from_slingshot_pseudotime_091202023/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_fitGAM_sce_object_10kvargenes_w_assocRes_startRes_09122023.rds")

sigGeneStart <- c("SOX2", "POU5F1", "HES5", "DLX5", "DLX1", "ONECUT2", "MYT1L", "HES6")

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts, gene = "REST")
plotGeneCount(crv, counts, gene = "DLX1")

install.packages("msigdbr")
library(msigdbr)

## C5 category is according to gene ontology grouping: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
geneSets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
### filter background to only include genes that we assessed.
geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
geneSets <- geneSets[geneSets$gene_symbol %in% names(sce),]
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
stats <- assocRes$waldStat_1
names(stats) <- rownames(assocRes)
genes <- rownames(assocRes)
BiocManager::install("fgsea")
library(fgsea)
eaRes <- fgsea(pathways = m_list, stats = stats, minSize = 10)
