#SCINA celltype software

install.packages('SCINA')
library(SCINA)
library(Seurat)
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(neuron.full.seu, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0("~/Desktop/", 'fetal_counts.mtx'))
df <- as.data.frame(counts_matrix)
write.csv(df, file='~/Desktop/fetal_counts.csv', quote=F, row.names=TRUE)


# write gene names
gene=rownames(counts_matrix)
genes <- as.data.frame(gene)
d <- cbind(genes$gene,df)
View(df)
exp=df


# Or .csv examples
exp=read.csv('your/path/to/example_expmat.csv',row.names=1,stringsAsFactors = F)
signatures=preprocess.signatures('/Volumes/Linda_Mac_backup_2/Rachel_scRNAseq/cell typing input files/pangleo_and_cellMarkerDB_cell_markers_r_input_version3.csv')

signatures=preprocess.signatures(marker_cellmarker_and_pangloDB)

results = SCINA(exp, marker_cellmarker_and_pangloDB, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=FALSE, log_file='SCINA.log')
View(results$cell_labels)
View(results$probabilities)

day2 <- AddMetaData(
  object = gcdata,
  metadata = results$cell_labels,
  col.name = 'SCINA_celltypes_7'
)

DimPlot(day2, reduction = "umap", group.by = "SCINA_celltypes_7", pt.size = .25)

write.csv(g, file='~/Desktop/day2_FindAllMarkers_3clusters.csv', quote=F, row.names=F)

write.table(gc@meta.data,file='~/Desktop/day7_celltype_metadata.txt', sep = "\t", quote=F,row.names=F,col.names=TRUE)

save(gcdata, file = "~/Desktop/gcdata_week4_celltype_annotations_final.RData")   

exp=df     

h <- read.csv("~/Desktop/test.csv", stringsAsFactors = FALSE)
                                 
                                 