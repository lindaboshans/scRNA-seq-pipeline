
devtools::install_github("JiaLiVUMC/scMRMA")

library(scMRMA)
load(system.file("data", "MouseBrain.Rdata", package = "scMRMA"))
result <- scMRMA(input = gcdata,
                 species = "Hu",
                 db = "panglaodb",
                 p = 0.05,
                 normalizedData = F,
                 selfDB = NULL,
                 selfClusters = NULL,
                 k=20)

result1 <- scMRMA(input=day2, species="Hu",selfClusters=Idents(day2))

CellType <- selfDefinedDatabase(file = "~/Desktop/markerlist_scMRMA.txt")


# UMAP plot
day2[["scMRMA_3"]] <- result$multiR$annotationResult[colnames(day2),ncol(result$multiR$annotationResult)]
DimPlot(day2,reduction = "umap",group.by = "scMRMA_3",label = TRUE,repel = TRUE)

