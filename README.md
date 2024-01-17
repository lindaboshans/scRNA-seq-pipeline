# scRNA-seq analysis Pipeline starting from 10x data

This repository contains a series of R scripts and Jupyter Notebooks for analyzing single-cell RNA-sequencing data from 10x raw data. It includes the following processes: Cell filtering and processing, cell clustering, sample integration for batch correction, cell annotation, similarity analysis, pseudotime analysis, trajectory inference, RNA velocity, and identification of key transcription factor regulomes. 

This is the basic workflow, with each step in separate scripts:

1. [Seurat](Seurat_dualSmad_AD_from_10x_processing.R)  R script contain code for starting from 10x raw data, creating seurat object, cell and mito filtering, and pre-processing steps.
2. [Seurat_analysis](Seurat_newer.R) contains key features of seurat, such as label transfer methods, integration, cell cycle regression, etc.
3. [Integration](Integration/) contains integration methods outside of seurat, such as CSS. 
4. [Gruffi_ribo_removal](gruffi_ribo_removal.R) contains the code to use the Gruffi package in order to remove stressed cells based on expression of key stress genes.
5. [Mrtree](mrtree.R) was used to determine the optimal number of clusters per sample.
6. The contents of [cell annotation](Cell%20Annotation/) contain all the types of annotation packages tested, some of which are cell marker based, and others which compare to gene expression of similar cell transcriptomes.
7. [Cell Trajectory and pseudotime analysis](Cell%20Trajectory%20and%20Pseudotime%20Analyses/) to determine the differentiation timeline of our cells.
8. [RNA velocity](Cell%20Trajectory%20and%20Pseudotime%20Analyses/all%20timepoints%20post%20CC%20regressed%20CSS%20integration%20RNA%20velocity%20scvelo.ipynb) 
9. [SCENIC](SCENIC.R) to identify key transcription factors involved in the differentiation process and their target genes.
10. [Brain Region Similarity Analyses](Fetal%20Brain%20Region%20Similarity/) to compare our cells to age-matched fetal cells during development. 
