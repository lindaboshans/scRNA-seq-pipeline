
#scCATCH
install.packages("scCATCH")
#or
# install devtools and install
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')

library(scCATCH)

load(paste0(system.file(package = "scCATCH"), "/extdata/mouse_kidney_203.rda"))
demo_geneinfo()

mouse_kidney_203 <- rev_gene(data = mouse_kidney_203, data_type = "data", species = "Mouse", geneinfo = geneinfo)
gc <- rev_gene(data = day2@assays$RNA@data, data_type = "data", species = "Human", geneinfo = geneinfo)

#character of seurat clusters
gc_cluster <- as.character(day2@meta.data$seurat_clusters)
#create scCATCH object
obj <- createscCATCH(data = gc, cluster = gc_cluster)

# find highly expressed genes
obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = c("Embryo", "Brain"), cell_min_pct = 0.1, logfc = 0.1)

#Evidence-based score and annotation for each cluster with findcelltype()

obj2 <- findcelltype(object = obj)

cellmatch_new <- marker_cellmarker_and_pangloDB
cellmatch_new <- cellmatch[cellmatch$species == "Human" & cellmatch$tissue %in% c( "Embryo", "Brain", "Fetal brain", "Embryonic brain", "Embryonic prefontal cortex", "Hippocampus", "Midbrain", "Dorsolateral prefrontal cortex", "Embryonic stem cell"), ]
obj1 <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj2 <- findcelltype(obj1)

r
