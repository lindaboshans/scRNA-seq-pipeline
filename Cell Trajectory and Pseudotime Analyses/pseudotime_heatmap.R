
library(data.table)
df <- fread("~/Desktop/alltimepoints_top2k_var_TFs_20231009.csv", header = FALSE)
Gene <- df$V1
Gene <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@assays$RNA@var.features
head(Gene)
meta <- fread("~/Desktop/gcdata_CC_regressed_CSS_integrated_metadata.txt")
summary(meta$slingshot_pseudotime)
dim(meta)

write.table(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@assays$RNA@var.features, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_var_features_2000.txt", quote=F, row.names=F, sep = '\t')
alltimepoints_top2k_var_TFs_20231009

# define quantile
ptime <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$phate_slingshot_pseudotime1
steps <- as.numeric(quantile(ptime, prob=seq(0, 1, by = 0.05)))
steps[15] <- 0.71882


write.csv(ptime, file = "~/Desktop/test.csv", row.names = TRUE)

# read count matrix
temp <- readMM("~/Desktop/gcdata_CC_regressed_CSS_integrated_counts.mtx")
temp <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@assays$RNA@counts
genenames <- read.table("~/Desktop/Desktop - Linda’s iMac/gcdata_CC_regressed_CSS_integrated_gene_names.csv")
genenames <- genenames$V1
genenames <- rownames(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@assays$RNA@counts)
rownames(temp) <- genenames
Mx <- temp[Gene, ]
dim(Mx); Mx[1:5, 1:5]


res <- matrix(data = 0, ncol = 20, nrow = length(Gene))
for (i in 1:(length(steps)-1)){
  print(paste(steps[i], steps[i+1]))
  cellindex <- ptime >=steps[i] & ptime < steps[i+1]
  res[, i] <- rowMeans(Mx[, cellindex])
}


gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins <- ""
res <- matrix(data = 0, ncol = 20, nrow = length(Gene))
for (i in 1:(length(steps)-1)){
  print(paste(steps[i], steps[i+1]))
  cellindex <- ptime >=steps[i] & ptime < steps[i+1]
  res[, i] <- rowMeans(Mx[, cellindex])
  gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins[cellindex] <- i
}

table(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins)


rownames(res) <- Gene
res <- res[rowSums(res) > 0, ]
write.csv(res, file = "~/Desktop/psuedotime_values_20_bins_res_2000_var_genes_09052023_mito_10and15.csv")

library(pheatmap)
color = colorRampPalette(c("blue", "white", "red"))(100)
pheatmap(res[, 1:20], cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 4, col = color, clustering_distance_rows = "correlation", clustering_method = "average", show_colnames = T)
pheatmap(res[, 1:20], cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 4, col = color, clustering_distance_rows = "correlation", clustering_method = "average", treeheight_row = 300)


# Compute the dissimilarity matrix
# df = the standardized data
res.dist <- dist(res, method = "euclidean")

#h clust to dteermine clustering method
res.hc <- hclust(d = res.dist, method = "complete")

# Compute cophentic distance
res.coph <- cophenetic(res.hc)

# Correlation between cophenetic distance and
# the original distance
cor(res.dist, res.coph)


library(ggplot2)
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins <- factor(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$bins, levels= 1:20) 
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$timepoint <- factor(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$timepoint, levels= unique(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$timepoint) )
ggplot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@meta.data, aes(x=bins, group=timepoint, fill=timepoint)) +
  geom_bar()


#extract genes from clusters via dendrogram branching
out <- pheatmap(res, cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 1, col = color, clustering_distance_rows = "correlation", clustering_method = "average", treeheight_row = 300)
#Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
d <-rownames(res[out$tree_row[["order"]],])
#You can also cut the tree at a pre-defined tree height, and extract the gene-to-cluster assignments at that height:
plot(out$tree_row, cex = 0.3)
abline(h=1, col="red", lty=2, lwd=2)#Cut the row (gene) dendrogram at a Euclidean distance dis-similarity of 8
three_clusters <- sort(cutree(out$tree_row, h=1.4))
four_clusters <- sort(cutree(out$tree_row, h=1))
five_clusters <- sort(cutree(out$tree_row, h=1.3))
seven_clusters <- sort(cutree(out$tree_row, h=1.65))

write.table(four_clusters, file = "~/Desktop/four_clusters_sai_psueodtime_heatmap_genes_top_2k_var_features_09052023.csv", sep = ",")
least_mat <- read.csv("~/Desktop/Desktop - Linda’s iMac/tradeseq_sig_DE_startvsEnd_down_all.csv", sep = ",", header = TRUE)

gene_list = least_mat$V1
names(gene_list) = least_mat$V1
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)
gene_list = least_mat$V1
gene_list = least_mat$Gene


BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

geneList = least_mat[,2]
names(geneList) = as.character(least_mat[,1])
geneList = sort(geneList, decreasing = TRUE)


geneList = rownames(res)
names(geneList) = as.character(least_mat[,1])
geneList = sort(geneList, decreasing = TRUE)


gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             scoreType = "pos",
             verbose = TRUE, 
             OrgDb = 'org.Hs.eg.db', 
             pAdjustMethod = "none")

rdotplot(gse, showCategory = 10)

ans.go <- enrichGO(gene = gene_list, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   universe = geneUniverse,
                   readable=TRUE,
                   keyType = "SYMBOL",
                   pvalueCutoff = 0.05)

tab.go <- as.data.frame(ans.go)
tab.go<- subset(tab.go, Count>5)
tab.go[1:5, 1:6]

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


c2.all <- read.gmt("~/Desktop/Desktop - Linda’s iMac/c2.cp.v2023.1.Hs.symbols.gmt")
c5.go <- read.gmt("~/Desktop/Desktop - Linda’s iMac/c5.go.v2023.1.Hs.symbols.gmt")
c8 <- read.gmt("~/Desktop/Desktop - Linda’s iMac/c8.all.v2023.1.Hs.symbols.gmt")

least_mat <- read.csv("~/Desktop/cluster 3 genes pseudotime heatmap sai.csv", sep = ",", header = FALSE)
least_mat <- read.csv("~/Desktop/tradeseq_startvsend_sig_DE_down_top_878_nonzero_padj.csv", sep = ",", header = TRUE)
least_mat <- read.csv("~/Desktop/tradeseq_startvsend_sig_DE_up_top_1000_nonzero_padj.csv", sep = ",", header = TRUE)
gene_list = least_mat$V1


ans.tf <- enricher(gene_list, TERM2GENE=c2.all)
tab.tf <- as.data.frame(ans.tf)
tab.tf<- subset(tab.tf, Count>5)
tab.tf[1:8,1:8]
write.table(tab.tf, file='~/Desktop/tradeseq_startVsEnd_up_top1000_sig_padj_genes_C2_CP_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf, file='~/Desktop/tradeseq_startVsEnd_down_top878_nonzero_sig_padj_genes_C2_CP_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf, file='~/Desktop/tradeseq_startVsEnd_up_top1000_nonzero_sig_padj_genes_C2_CP_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf, file='~/Desktop/cluster3_slingshot_pseudotime_bins_sai_methodC2_CP_mito10and15_top_2k_var_genes_stress_09042023.csv', sep = ",")
dotplot(ans.tf, showCategory = 10)

ans.tf1 <- enricher(gene_list, TERM2GENE=c5.go)
tab.tf1 <- as.data.frame(ans.tf1)
tab.tf1<- subset(tab.tf1, Count>5)
tab.tf1[1:8,1:8]
write.table(tab.tf1, file='~/Desktop/tradeseq_startVsEnd_up_top1000_nonzero_sig_padj_genes_C5_GO_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf1, file='~/Desktop/tradeseq_startVsEnd_down_top878_nonzero_sig_padj_genes_C5_GO_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf1, file='~/Desktop/tradeseq_startVsEnd_up_top1000_nonzero_sig_padj_genes_C5_GO_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf1, file='~/Desktop/cluster3_slingshot_pseudotime_bins_sai_method_C5_GO_mito10and15_top_2k_var_genes_stress_09042023.csv', sep = ",")
dotplot(ans.tf1, showCategory = 10)

ans.tf2 <- enricher(gene_list, TERM2GENE=c8)
tab.tf2 <- as.data.frame(ans.tf2)
tab.tf2<- subset(tab.tf2, Count>5)
tab.tf2[1:8,1:8]
write.table(tab.tf2, file='~/Desktop/tradeseq_startVsEnd_up_top1000_sig_padj_genes_C8_celltype_signature_gene_sets_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf2, file='~/Desktop/tradeseq_startVsEnd_down_top878_sig_padj_genes_C8_celltype_signature_gene_sets_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf2, file='~/Desktop/tradeseq_startVsEnd_up_top1000_sig_padj_genes_C8_celltype_signature_gene_sets_mito10and15_05312023.csv', sep = ",")
write.table(tab.tf2, file='~/Desktop/cluster3_slingshot_pseudotime_bins_sai_method_C8_celltype_mito10and15_top_2k_var_genes_stress_09042023.csv', sep = ",")
dotplot(ans.tf2, showCategory = 10)


dotplot(ans.tf, showCategory = 10)



pheatmap(data[, 3:9], cluster_rows = T, cluster_cols = F, scale = "row", fontsize_row = 4, col = color, clustering_distance_rows = "correlation", clustering_method = "average")

