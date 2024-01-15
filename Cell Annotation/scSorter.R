

#scSorter for celltyping clusters

install.packages("remotes")
remotes::install_github("hyguo2/scSorter")
library(scSorter)

library(Seurat)
#expr_obj = CreateSeuratObject(expr)
#or
load("~/Desktop/gcdata_day7_celltype_annotations_final.RData")


#marker annotation file should look like below 
head(anno)
#>                 Type Marker Weight
#> 1  Endothelial_cells Pecam1      2
#> 2  Endothelial_cells   Cdh5      2
#> 3  Endothelial_cells    Kdr      2
#> 4             Immune  Ptprc      2
#> 5 Pancreatic_A_cells    Gcg      2
#> 6 Pancreatic_A_cells   Mafb      2



expr_obj <- NormalizeData(gc4, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)

#next, we choose highly variable genes by FindVariableFeatures() function. 
#We also filter out genes with non-zero expression in less than 10% of total cells.

expr_obj <- FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
topgenes <- head(VariableFeatures(expr_obj), 2000)

expr = GetAssayData(expr_obj)
topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

#last, we subset the preprocessed expression data. Now, we are ready to run scSorter.

picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]



