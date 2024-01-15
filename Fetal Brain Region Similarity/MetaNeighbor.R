
library(MetaNeighbor)
library(SummarizedExperiment)
source("2016-11-03-runMetaNeighbor.R")
#load("MetaNeighbor_sample_data.Rdata")

fetal <- subset(wk4_fetal_brain_seurat_integrated, subset = structure_collapsed == "claustrum", invert = TRUE)


sce <- as.SingleCellExperiment(fetal)
dim(sce)
head(colData(sce))
table(sce$"ident", sce$structure_collapsed)

global_hvgs = variableGenes(dat = sce, exp_labels = sce$dataset, min_recurrence=1)
#global_hvgs = wk4_fetal_brain_seurat_integrated@assays$integrated@var.features
length(global_hvgs)

#ideal var genes is 600 to 1000. 
global_hvgs <- global_hvgs[1:1000]
length(global_hvgs)


aurocs = MetaNeighborUS(var_genes = global_hvgs,
                        dat = sce,
                        study_id = sce$dataset,
                        cell_type = sce$structure_collapsed,
                        fast_version = TRUE)

plotHeatmap(aurocs, cex = .5)

topHits(aurocs, dat = sce, study_id = sce$dataset,
        cell_type = sce$structure_collapsed, threshold = 0.9)

best_hits = MetaNeighborUS(var_genes = global_hvgs,
                           dat = sce,
                           study_id = sce$dataset,
                           cell_type = sce$structure_collapsed,
                           fast_version = TRUE,
                           one_vs_best = TRUE, symmetric_output = FALSE)
plotHeatmap(best_hits, cex = 0.5)




