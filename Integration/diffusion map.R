
#make single cell object from seurat
library(Seurat) 
library(SingleCellExperiment)
load("~/Desktop/gdata_CC_regressed_CSS_integrated_1000_union_genes_res_0.6_wslingshot_pseudotime.RData")

sce <- as.SingleCellExperiment(gcdata_CC_regressed_CSS_integrated)
#this has the cell classification
table(sce$ident)

#this step may take a long time (days) or not finish. It is recommend to send it to the cluster as a script that reads the Seurat or the single cell object, runs DiffusionMap, and saves the object. 
library(destiny)
dm <- DiffusionMap(sce, verbose = TRUE)

#  Prepare a counts matrix with labeled rows and columns. 
deng <- logcounts(sce)  # access log-transformed counts matrix
cellLabels <- sce$cell.name
colnames(deng) <- cellLabels
mtx <- as.matrix(deng)

# Make a diffusion map.
dm <- DiffusionMap(t(mtx))

# Optional: Try different sigma values when making diffusion map.
# dm <- DiffusionMap(t(deng), sigma = "local")  # use local option to set sigma
# sigmas <- find_sigmas(t(deng), verbose = FALSE)  # find optimal sigma
# dm <- DiffusionMap(t(deng), sigma = optimal_sigma(sigmas))  

# Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2). 
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()