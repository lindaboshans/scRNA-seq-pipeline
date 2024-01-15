#scRNA-seq analysis Pipeline starting from 10x data

This is the basic workflow, with each step in separate scripts:

1. Seurat_dualSmad_AD_from_10x_processing.R R script contain code for starting from 10x raw data, creating seurat object, cell and mito filtering, and pre-processing steps.
2. Suerat_newer contains key features of seurat, such as label transfer methods, integration, cell cycle regression, etc.
3. Gruffi_ribo_removal contain the code to use the Gruffi package in order to remove stressed cells based on expression of key stress genes.
4. Mrtree was used to determine the optimal number of clusters per sample.
5. The contents of cell annotation contain all the types of annotation packages tested, some of which are cell marker based, and others which compare to gene expression of similar cell transcriptomes.
6. Trajectory inference and pseudotime analysis to determien the differentiation timeline of our cells.
7. 
