install.packages("devtools")
library(devtools)

# installing/loading the package:
if(!require(installr)) {
  install.packages("installr"); 
  require(installr)
} #load / install+load installr

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "seqLogo"))
install.packages("phangorn")
library(phangorn)
remotes::install_github("KlausVigo/phangorn")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

‘Biobase’, ‘IRanges’, ‘S4Vectors’, ‘SummarizedExperiment’ 

devtools::install_github("YosefLab/SymSim")
library(SymSim)

BiocManager::install("ggtree")
BiocManager::install("SingleCellExperiment")
library(ggtree)
library(SingleCellExperiment)

devtools::install_github("pengminshi/mrtree", force = TRUE)

library(mrtree)


set.seed(42)


week4_slc32a1_only <- subset(gcdata, subset = SLC32A1 >0)

metadata = wk4_slc32a1_only@meta.data
ref.labels = wk4_slc32a1_only@meta.data$seurat_clusters
counts_matrix <- GetAssayData(wk4_slc32a1_only, assay='RNA', slot='counts')
df <- as.matrix(counts_matrix)

metadata = GW20_brain_regions_combined@meta.data
ref.labels = GW20_brain_regions_combined@meta.data$region.ident
counts_matrix <- GetAssayData(GW20_brain_regions_combined, assay='RNA', slot='counts')
df <- as.matrix(counts_matrix)


metadata = wk4_interneuron@meta.data
ref.labels = wk4_interneuron@meta.data$seurat_clusters
counts_matrix <- GetAssayData(wk4_interneuron, assay='RNA', slot ='counts')
df <- as.matrix(counts_matrix)


d# specify the resolution parameters
resolutions = seq(0.1, sqrt(3), 0.1)^2

# alternatively and preferrably, we provide a sampling tool to 
# sample resolution parameters to uniformly cover different scales
A = seurat_get_nn_graph(counts=wk4_interneuron@assays$RNA@counts, metadata=metadata, npc=20)
resolutions = modularity_event_sampling(A=A, n.res=30, 
                                        gamma.min=0.01, gamma.max=2.5) # sample based on the similarity matrix

# clustering using Suerat 
gcdata.out = sc_clustering.seurat(counts=df, resolutions=resolutions, 
                                  metadata=metadata, npcs=20,
                                  return.seurat.object=TRUE, vars.to.regress=NULL,
                                  find.variable.features=FALSE,verbose=FALSE)

# initial cluster tree from Seurat flat clustering
plot_clustree(labelmat=gcdata.out$seurat.clusters,  prefix ='RNA_snn_res.', suffix = '0.4',
              ref.labels = ref.labels, plot.ref = FALSE)
#> Warning: The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
#> Please use the `.add` argument instead.

out = mrtree(gcdata.out$obj, n.cores = 4, consensus=FALSE, augment.path=FALSE, max.k = 20)
# if there are few partitions per k, within resolution consensus step can speed up the algorithm
# weight per sample is encoraged if the classes are imbalanced

plot_tree(labelmat=out$labelmat.mrtree, ref.labels=ref.labels, plot.piechart = TRUE,
          node.size = 0.1, tip.label.dist = 10, bottom.margin=20 )

save(out, file = "~/Desktop/day7_mrtree_results.RData")

ks.flat = apply(out$labelmat.flat, 2, FUN=function(x) length(unique(x)))
ks.mrtree = apply(out$labelmat.mrtree, 2, FUN=function(x) length(unique(x)))
amri.flat = sapply(1:ncol(out$labelmat.flat), function(i) AMRI(out$labelmat.flat[,i], ref.labels)$amri)
amri.flat = aggregate(amri.flat, by=list(k=ks.flat), FUN=mean)
amri.recon = sapply(1:ncol(out$labelmat.mrtree), function(i) AMRI(out$labelmat.mrtree[,i], ref.labels)$amri)

df = rbind(data.frame(k=amri.flat$k, amri=amri.flat$x, method='Seurat flat'), 
           data.frame(k=ks.mrtree, amri=amri.recon, method='MRtree'))
ggplot2::ggplot(data=df, aes(x=k, y=amri, color=method)) + geom_line() + theme_bw()
#> Warning: Removed 2 row(s) containing missing values (geom_path).
#> 
stab.out = stability_plot(out)
stab.out$plot


week4_slc32a1_only

gcdata_mito_filtered <- RunUMAP(gcdata_mito_filtered, assay = "RNA", dims = 1:20)
gcdata_mito_filtered <- FindNeighbors(gcdata_mito_filtered)
gcdata_mito_filtered <- FindClusters(gcdata_mito_filtered, resolution = 0.15)

# Visualization
DimPlot(gcdata_mito_filtered, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(gcdata_CC_regressed_CSS_integrated_slc32a1_and_gad1_and_gad2_week4, reduction = "umap", group.by = "structure")

p1 + p2
p1




