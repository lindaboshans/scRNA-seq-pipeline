
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("kstreet13/slingshot")

BiocManager::install("scran")
BiocManager::install("GenomeInfoDb")
BiocManager::install("DelayedMatrixStats")
library(DelayedMatrixStats)
library(slingshot)
library(scran)

suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(igraph)
})

load("~/Desktop/gcdata_with_dualSMAD_mito_10_and_15_and20_CC_regressed_CSS_integrated_1000_union_genes_res_0.9_w_slingshot_pseudotime_20230601.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Rachel_scRNAseq/RData_files/gcdata_LB.RData")

load("~/Desktop/gcdata_CC_regressed_CSS_integrated.RData")
save(gcdata_CC_regressed_CSS_integrated, file = "~/Desktop/gcdata_with_dualSMAD_mito_10_and_15_and20_CC_regressed_CSS_integrated_1000_union_genes_res_0.9_w_slingshot_pseudotime_20230601.RData")

load("~/Desktop/gcdata_percent_mito_10_and_15_CC_regressed_CSS_integrated_union_top1000genes_res_0.9_w_slingshot_pseudotime_20230518.RData")
load("~/Desktop/gcdata_with_dualSMAD_mito_10_and_15_and20_CC_regressed_CSS_integrated_1000_union_genes_res_0.9_20230601.RData")


# Assuming your Seurat object is named "seuratObj"
# Assuming the column you want to change from character to factor is named "timepoint"

# Convert the column from character to factor
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated$timepoint <- as.factor(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated$timepoint)
gcdata_CC_regressed_CSS_integrated$seurat_clusters <- as.factor(gcdata_CC_regressed_CSS_integrated$seurat_clusters)



# Save the objects as separate matrices for input in slingshot
#dimred <- gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated@reductions$scVI_MDE@cell.embeddings
#clustering <- gcdata_CC_regressed_CSS_integrated$RNA_snn_res.1
#counts <- as.matrix(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated@assays$RNA@counts[gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated@assays$RNA@var.features, ])

#counts <- as.matrix(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated@assays$RNA@counts)



#sds <- slingshot(Embeddings(gcdata_CC_regressed_CSS_integrated, "scVI_MDE"), clusterLabels = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_CSS_integrated$timepoint, 
                 #start.clus = "48h_pD", stretch = 0)


#cell_colors <- cell_pal(gcdata_CC_regressed_CSS_integrated$timepoint, brewer_pal("qual", "Set2"))
#cell_colors_clust <- cell_pal(gcdata_CC_regressed_CSS_integrated$seurat_clusters, hue_pal())



#start here

sce <- as.SingleCellExperiment(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed,  assay = "RNA")
sce <- slingshot(sce, clusterLabels = 'timepoint', reducedDim = 'PHATE')
summary(sce$slingPseudotime_1)

#plot(reducedDims(sce)$UMAP_CSS, col = col[cut(sce$slingPseudotime_1,breaks=5)], pch=16, asp = 1, cex = 0.2)
#lines(SlingshotDataSet(sce), lwd=2)

#plotGenePseudotime(sce, 'NR2F1')


library(grDevices)
install.packages('Polychrome')
library(Polychrome)
library(RColorBrewer)
#colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
#my_color <- createPalette(10, c("#010101", "#ff0000"), M=1000)
#col <- brewer.pal(8,'Set1')
#col <- brewer.pal(6,'Dark2')

col <-  c(
  RColorBrewer::brewer.pal(5,'Set1')[1],
  RColorBrewer::brewer.pal(5,'Set1')[5],
  RColorBrewer::brewer.pal(6,'Dark2')[6],
  RColorBrewer::brewer.pal(5,'Set1')[3],
  RColorBrewer::brewer.pal(5,'Set1')[2],
  RColorBrewer::brewer.pal(5,'Set1')[4]
)

#RColorBrewer::brewer.pal(7,'Set1')[7]

my_cols <- c('48h_pD'='#E41A1C','7d_postDox'='#FF7F00','1wk_glia'='#E6AB02','2wk_glia'='#4DAF4A','3wk_glia'='#377EB8',
             '4wk_glia'='#984EA3','dualSMAD'='#A65628')

my_cols <- c('48h_pD'='#E41A1C','7d_postDox'='#FF7F00','1wk_glia'='#E6AB02','2wk_glia'='#4DAF4A','3wk_glia'='#377EB8',
             '4wk_glia'='#984EA3')

#clus <- c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia', 'dual_SMAD')
#colvec <- col[factor(clus, levels = c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia', 'dualSMAD'))]


install.packages("viridis")

#colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#plotcol <- colors[cut(sce$slingshot_pseudotime, breaks=100)]
#plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1, cex = .3, main = 'Lineage 1') 
#lines(SlingshotDataSet(sce), lwd=2, col='black')


#par(mfrow = c(1, 2))
#colfunc <- colorRampPalette(colors)
#legend_image <- as.raster(matrix(rev(colfunc(20)), ncol=1))
#plot(reducedDims(sce)$UMAP_CSS, col = plotcol, pch=16, asp = 1, cex = .25, main = 'slingshot Pseudotime')
#lines(SlingshotDataSet(sce), lwd=2, col='black')
#text(x=-10, y = seq(10.6, 6.4,l=5), labels = rev(seq(0,30,l=5)), cex = .6)
#mtext('Pseudotime', -1.5, -.125, cex=.6, adj=.1, font=2)  
#rasterImage(legend_image, -11.6, 6.3, -10.5, 10.6)


library(viridis)

#plot(reducedDims(sce)$UMAP, col = col[as.character(sce$slingPseudotime_1)], pch=16, asp = 1, cex = .25)
#lines(SlingshotDataSet(sce), lwd=2, col='black')


#plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$slingPseudotime_1], pch=16, asp = 1)
#lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')


slingLineages(sce)
summary(sce$slingPseudotime_1)
par(mfrow = c(1, 2))
plot(reducedDims(sce)$PHATE, col = col[sce$timepoint], pch=16, asp = 1, cex = .25)
legend('topright', title = 'timepoint', col = col, legend=c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia', 'dual_SMAD'), pch=16, cex = 0.5)
lines(SlingshotDataSet(sce), lwd=1.5, type = 'curves', col = 'black')
lines(SlingshotDataSet(sce), lwd=1.5, type = 'lineages', col = 'black')

#is.liblegend("topleft", legend = my_color[levels(sce$timepoint)],  
       #fill = my_color[levels(sce$timepoint)], inset=0.05)
#legend("topleft", legend = my_color[sce$timepoint],  
       #fill = my_color[sce$timepoint], inset=0.05)
#legend("topright", legend=unique(sce$timepoint), title="Group", col = brewer.pal(9,'Set1'))
#sds <- as.SlingshotDataSet(sce)
#plot(rd, col = 'grey50', asp = 1)
#lines(curve1, lwd=2)


#lin1 <- getLineages(sce, start.clus = '1')
#lin1

#plot(sce, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
#lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

## Plotting the pseudotime inferred by slingshot by cell types

library(ggbeeswarm)



pseudo_1 <- sce$slingPseudotime_1
pseudo_2 <- sce$slingPseudotime_2
timepoint <- as.data.frame(sce$timepoint)
umap <- reducedDims(sce)$PHATE
UMAP <- as.data.frame(umap)
df <- colData(sce)

df <- cbind(pseudo_1, timepoint)
df1 <- cbind(pseudo_2, df)
df2 <- cbind(umap, df1)
df2 <- cbind(umap, df)


#scale values to scale of 0 to 1
normalized = (x-min(x))/(max(x)-min(x))

#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#apply Min-Max normalization to first four columns in iris dataset
new_pseudo1 <- min_max_norm(sce$slingPseudotime_1)
sce$slingPseudotime_1 <- new_pseudo1

norm = (.06-0)/(.08639-0)

gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed <- AddMetaData(object = gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, metadata = new_pseudo1, col.name = 'phate_slingshot_pseudotime1')
gcdata_CC_regressed_CSS_integrated<- AddMetaData(object = gcdata_CC_regressed_CSS_integrated, metadata = pseudo_2, col.name = 'slingshot_pseudotime2')

ggplot(df2, aes(df2$pseudo_1, df2$`sce$timepoint`)) + geom_quasirandom(groupOnX = FALSE) +
  geom_point(aes(color = colors),
             alpha = 0.5) +
  theme_minimal() + labs(colour = "Pseudotime")


#this one works for pseudotime
ggplot(UMAP, aes(x = new_pseudo1, y = sce$timepoint, 
                                             colour = sce$timepoint)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_colour_manual(values = col) + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime1")


slingPseudotime(sce)

ggplot(df2, aes(x = sce$slingPseudotime_1, y = sce$timepoint, 
                                             colour = sce$timepoint)) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime_1")


embedded <- slingCurves(embedCurves(sce, "UMAP"))

gg <- plotUMAP(sce,text_by = "timepoint", 
               colour_by = "timepoint", point_size = 0.5)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=Dim.1, y=Dim.2), size=1.2)
}

# Pseudotime densities (by spatial)

new_pseudo1 <- as.data.frame(new_pseudo1)
row.names(new_pseudo1) <- colnames(sce)

x <- na.omit(slingPseudotime(sce)[colData(sce)$timepoint == "48h_pD", 1])

ds <- list(day2 = density(slingPseudotime(sce)[colData(sce)$timepoint == "48h_pD", 1], na.rm = TRUE),
            wk2 = density(slingPseudotime(sce)[colData(sce)$timepoint == "2wk_glia", 1], na.rm = TRUE),
            wk1 = density(slingPseudotime(sce)[colData(sce)$timepoint == "1wk_glia", 1], na.rm = TRUE),
            day7 = density(slingPseudotime(sce)[colData(sce)$timepoint == "7d_postDox", 1], na.rm = TRUE),
            wk4 = density(slingPseudotime(sce)[colData(sce)$timepoint == "4wk_glia", 1], na.rm = TRUE),
            wk3 = density(slingPseudotime(sce)[colData(sce)$timepoint == "3wk_glia", 1], na.rm = TRUE))
          

ds <- list(day2 = density(new_pseudo1[colData(sce)$timepoint == "48h_pD", 1], na.rm = TRUE),
           wk2 = density(new_pseudo1[colData(sce)$timepoint == "2wk_glia", 1], na.rm = TRUE),
           wk1 = density(new_pseudo1[colData(sce)$timepoint == "1wk_glia", 1], na.rm = TRUE),
           day7 = density(new_pseudo1[colData(sce)$timepoint == "7d_postDox", 1], na.rm = TRUE),
           wk4 = density(new_pseudo1[colData(sce)$timepoint == "4wk_glia", 1], na.rm = TRUE),
           wk3 = density(new_pseudo1[colData(sce)$timepoint == "3wk_glia", 1], na.rm = TRUE))



wk4 = density(slingPseudotime(sce)[colData(sce)$timepoint == "4wk_glia", 1], na.rm = TRUE),)

col <-  c(
  RColorBrewer::brewer.pal(5,'Set1')[1],
  RColorBrewer::brewer.pal(5,'Set1')[5],
  RColorBrewer::brewer.pal(6,'Dark2')[6],
  RColorBrewer::brewer.pal(5,'Set1')[3],
  RColorBrewer::brewer.pal(5,'Set1')[2],
  RColorBrewer::brewer.pal(5,'Set1')[4]
)

xlim <- range(c(ds$day2$x, ds$wk4$x))
ylim <- range(c(ds$day2$y, ds$wk3$y))
plot(xlim, ylim, col = "white", xlab = "Pseudotime", ylab = "Density")
polygon(c(min(ds$day2$x),  ds$day2$x, max(ds$day2$x)), c(0, ds$day2$y, 0),
        col = alpha(brewer.pal(6, "Set1")[1], alpha = .5))
polygon(c(min(ds$day7$x), ds$day7$x, max(ds$day2$x)), c(0, ds$day7$y, 0),
        col = alpha(brewer.pal(6, "Set1")[5], alpha = .5))
polygon(c(min(ds$wk1$x), ds$wk1$x, max(ds$day2$x)), c(0, ds$wk1$y, 0),
        col = alpha(brewer.pal(6, "Dark2")[6], alpha = .5))
polygon(c(min(ds$wk2$x), ds$wk2$x, max(ds$day2$x)), c(0, ds$wk2$y, 0),
        col = alpha(brewer.pal(6, "Set1")[3], alpha = .5))
polygon(c(min(ds$wk3$x), ds$wk3$x, max(ds$day2$x)), c(0, ds$wk3$y, 0),
        col = alpha(brewer.pal(6, "Set1")[2], alpha = .5))
polygon(c(min(ds$wk4$x), ds$wk4$x, max(ds$day2$x)), c(0, ds$wk4$y, 0),
        col = alpha(brewer.pal(6, "Set1")[4], alpha = .5))
polygon(c(min(ds$dualSMAD$x), ds$dualSMAD$x, max(ds$dualSMAD$x)), c(0, ds$dualSMAD$y, 0),
        col = alpha(brewer.pal(7, "Set1")[7], alpha = .5))

legend("topleft",legend=c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia'), 
       fill = alpha(col, alpha = .5), bty = "n")

plot(reducedDims(sce)$UMAP_CSS[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(3, 4)[factor(colData(sce)$timepoint)][shuffle], alpha = .5))
lines(SlingshotDataSet(sce), type = 'lineages', show.constraints = TRUE)
legend("topright", pch = 16, col = 3:4, bty = "n", legend = levels(factor(colData(sce)$timepoint)))
layout(1)
par(mar = c(5, 4, 4, 2) + .1)


#projecting cells onto existing trajectories
# our original PseudotimeOrdering
pto <- sce$slingshot

# simulate new cells in PCA space
newPCA <- reducedDim(sce, 'PHATE') + rnorm(2*ncol(sce), sd = 2)
embedded <- embedCurves(sce, newPCA)

# project onto trajectory
newPTO <- slingshot::predict(pto, newPCA)

color_names <- colors()  # Correctly retrieve the list of color names
# Ensure that the number of breaks in cut does not exceed the length of color_names
newplotcol <- color_names[cut(slingPseudotime(newPTO)[, 1], breaks=min(100, length(color_names)))]


newplotcol <- colors[cut(slingPseudotime(newPTO)[,1], breaks=100)]
plot(reducedDims(sce)$PCA, col = 'grey', bg = 'grey', pch=16, asp = 1,
     xlim = range(newPCA[,1]), ylim = range(newPCA[,2]))
lines(SlingshotDataSet(sce), lwd=2, col = 'black')
points(slingReducedDim(newPTO), col = col[sce$timepoint], pch = 16)
  
par(mfrow = c(1, 2))
plot(reducedDims(sce)$PHATE, col = col[sce$timepoint], pch=16, asp = 1, cex = .25)
legend('topright', title = 'timepoint', col = col, legend=c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia', 'dual_SMAD'), pch=16, cex = 0.5)
lines(SlingshotDataSet(sce), lwd=2, type = 'curves', col = 'black')
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
  

#project pseudotime values onto different reducedDIM
plot(umap, asp=1, col = 'grey75', pch = 16, cex = .5)
points(umap, col = hcl.colors(100)[cut(slingPseudotime(sce)[,1], 100)], pch = 16, cex = .5)
legend('bottomright', title = 'Pseudotime', col = hcl.colors(5), legend=c('0','.25', '.5', '.75', '1'), pch=16, cex = .7)


  #Keep only cells with UMAP coordinates (using CATALYST function)
  sce_sling <- filterSCE(sce, complete.cases(reducedDim(sce))) #remove NA
  
  #Prepare clusterLabels
  clusters <- sce$cluster_id
  levels(clusters) <- cluster_codes(sce)$cluster_annotation
  
  #Run slingshot
  sce_sling <- slingshot(sce_sling, clusterLabels = clusters, 
                         start.clus = "C2", stretch = 0)
  
  #Create dataframe of pseudotimes
  pt <- setNames(as.data.frame(slingPseudotime(sce, na = FALSE)), c("slingPseudotime_1", "slingPseudotime_2"))
  
  #Run slingCurves
  curve2 <- slingCurves(sce)[[2]]
  cv1 <- setNames(curve1$s[curve1$ord,] %>% as.data.frame(), c("UMACSS_1", "UMACSS_2"))
  
  #Create dataframe with UMAP coordinates
  umap_df <- setNames(as.data.frame(reducedDim(sce)), c("UMACSS_1", "UMACSS_2"))
  
  #Create dataframe with UMAP coordinates and pseudotimes
  df <- cbind(umap_df, pt)
  
  #Plot according to pseudotime values
  p1 <- ggplot(df2, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = df2$pseudo_1),
               alpha = 0.5) +
    scale_colour_viridis_c() +
    theme_minimal() + labs(colour = "Pseudotime")
  
  tiff("./UMAPpseudo.tiff", width = 5*900, height = 5*900, res = 300, pointsize = 5)     
  p1
  dev.off()
  
  #Extract UMAP_1 and UMAP_2 and plot 
  UMAP_1 <- p1$data$UMACSS_1
  UMAP_2 <- p1$data$UMACSS_2
  
  tiff("./sling.tiff", width = 5*400, height = 5*300, res = 300, pointsize = 5)     
  p1 + geom_path(aes(x = UMAP_1, y = UMAP_2), data = cv1,
                 col = "black", linewidth = 1, arrow = arrow(), lineend = "round") 
  dev.off()
  
  
  
  #Which cells are in which lineage? Here we plot the pseudotime values for each lineage.
  install.packages("viridis")
  library(viridis)
  
  
  nc <- 2
  pt <- slingPseudotime(sce)
  nms <- colnames(pt)
  nr <- ceiling(length(nms)/nc)
  pal <- viridis(100, end = 0.95)
  par(mfrow = c(nr, nc))
  for (i in nms) {
    colors <- pal[cut(pt[,i], breaks = 100)]
    plot(reducedDims(sce)$UMAP, col = colors, pch = 16, cex = 0.2, main = i)
    legend('bottomleft', title = 'pseudotime', col = c("#440154FF", "#365C8DFF", "#22A785FF", "#DEE318FF" ), legend=c('0','10','20', '30'), pch=16)
    lines(SlingshotDataSet(sce), lwd = 2, col = 'black', type = 'curves')
  }
  
  colors <- pal[cut(pt[,"Lineage1"], breaks = 100)]
  plot(reducedDims(sce)$UMAP_CSS, col = colors, pch = 16, cex = 0.2, main = "Lineage1")
  legend('bottomright', title = 'pseudotime', col = c("#440154FF", "#365C8DFF", "#22A785FF", "#DEE318FF" ), legend=c('0','10','20', '30'), pch=16)
  lines(SlingshotDataSet(sce), lwd = 2, col = 'black', type = 'curves', linInd = 1)
  
  colors <- pal[cut(pt[, "Lineage2"], breaks = 100)]
  plot(reducedDims(sce)$UMAP_CSS, col = colors, pch = 16, cex = 0.2, main = "Lineage2")
  legend('bottomright', title = 'pseudotime', col = c("#440154FF", "#365C8DFF", "#22A785FF", "#DEE318FF" ), legend=c('0','10','20', '30'), pch=16)
  lines(SlingshotDataSet(sce), lwd = 2, col = 'black', type = 'curves', linInd = 2)


  plot(reducedDims(sce)$UMAP_CSS, col = colors[sce$timepoint], pch=16, asp = 1, cex = .2)
  legend('bottomright', title = 'timepoint', col = color(6), legend=c('48h_pD','7d_postDox','1wk_glia', '2wk_glia','3wk_glia', '4wk_glia'), pch=16)
  legend('bottomright', title = 'timepoint', col = col, legend=c(sce$timepoint), pch=16)
  lines(SlingshotDataSet(sce), lwd=2, type = 'curves', col = 'black')
  lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')


  install.packages("gam")
  library(gam)
  # Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
  # Identify the variable genes by ranking all genes by their variance.
  
  Y <- log2(counts(sce) + 1)
  var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
  Y <- Y[var1K, ]  # only counts for variable genes
  
  # Fit GAM for each gene using pseudotime as independent variable.
  t <- sce$slingPseudotime_1
  #because we scaled values, using seurat object to pull values
  t <- as.data.frame(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$phate_slingshot_pseudotime1)
  t <- t$`gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed$phate_slingshot_pseudotime1`
  
  
  gam.pval <- apply(Y, 1, function(z){
    d <- data.frame(z=z, t=t)
    tmp <- gam(z ~ lo(t), data=d)
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  
  write.csv(gam.pval, file = "~/Desktop/gam_pval.csv", row.names =TRUE)
  
  # Identify genes with the most significant time-dependent model fit.
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]  
  
  # Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
  BiocManager::install("clusterExperiment")
  
  require(clusterExperiment)
  heatdata <- as.matrix(gcdata_CC_regressed_CSS_integrated@assays$RNA@counts[rownames(gcdata_CC_regressed_CSS_integrated@assays$RNA@counts) %in% topgenes, order(t, na.last = NA)])
  write.csv(heatdata, file = "~/Desktop/heatdata.csv", row.names =F)
  heatclus <- gcdata_CC_regressed_CSS_integrated$timepoint[order(t, na.last = NA)]
  png(paste0("~/Desktop/", "heatmap_pseudotime_genes_RNA_counts_TFs2.png"), width=10, height=10, units = "in", res=200)
  ce <- ClusterExperiment(heatdata, heatclus, transform = log1p)
  clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', cexRow = 1.5, fontsize = 15, colorScale = seqPal2)
  dev.off()

  topgenes1 <- read.delim("~/Desktop/top_tfs.csv", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  topgenes1 <- names(sort(topgenes, decreasing = FALSE))

  #add quations to list 
  topgenes <- gsub("\b", "'", topgenes1)

topgenes <- c("HMGA1", "FOS", "TSC22D1", "BARX1", "POU5F1", "SOX4", "JUND", "YBX1", "JUN", "ZNF428", "HES6", "TERF1", "GTF2I", "LIN28A", "SON", "PA2G4", "BPTF", "HMGN3", "ZNF292", "SOX11", "KLF6", "ATF4", "DRAP1", "BAZ2B", "XBP1", "PIN1", "UNCX", "KDM5B", "CAMTA1", "CXXC5", "CSRNP3", "DLX5", "YY1", "PBX1", "TFDP2", "FOXP1")

#identify DE genes alon trajectory
pseudo <- testPseudotime(sce, pseudotime=sce$slingshot[,1])[[1]]
pseudo$SYMBOL <- rowData(sce)$SYMBOL
pseudo_ordered <- pseudo[order(pseudo$p.value),]

# Making a copy of our SCE and including the pseudotimes in the colData.
sce2 <- sce
sce2$TSCAN.first <- pathStat(sce$slingshot)[,1]
#sce2$TSCAN.second <- pathStat(sce$slingshot)[,2]

# Discarding the offending cluster.
#discard <- "7"
#keep <- colLabels(sce.nest)!=discard
#sce.nest2 <- sce.nest2[,keep]

# Testing against the first path again.
pseudo <- testPseudotime(sce2, pseudotime=sce2$TSCAN.first)
pseudo$SYMBOL <- rowData(sce.nest2)$SYMBOL
sorted <- pseudo[order(pseudo$p.value),]

up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)

best <- head(up.left, 10)
plotExpression(sce2, features=best, swap_rownames="SYMBOL",
               x="TSCAN.first", colour_by="label")

up.right <- sorted[sorted$logFC > 0,]
head(up.right, 10)

#heatmap
on.first.path <- !is.na(sce2$TSCAN.first)
plotHeatmap(sce2[,on.first.path], order_columns_by="TSCAN.first", 
            colour_columns_by="label", features=head(up.right, 50),
            center=TRUE, swap_rownames="SYMBOL")
