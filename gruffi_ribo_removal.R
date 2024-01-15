BiocManager::install("DOSE")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("sparseMatrixStats")
BiocManager::install("biomaRt")

devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", ref = 'v2.3.4', upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
BiocManager::install("EnhancedVolcano")
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)

devtools::install_github(repo = "jn-goe/gruffi", upgrade = T)

options(rgl.useNULL = TRUE)
library(rgl)
library(gruffi)
library(ggExpress)
library(biomaRt)


if(is.null(sum(grepl(".reassigned", Seurat.utils::GetClusteringRuns(AD3))))) 
  Seurat.utils::GetClusteringRuns(AD3)[1] 
    else Seurat.utils::GetClusteringRuns(AD3)[grepl(".reassigned", Seurat.utils::GetClusteringRuns(AD3))]

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering
go4 <- "GO:0022008" # Neurogenesis, negative filtering

mature_INs_slc32a1_or_gad_only <- aut.res.clustering(obj = mature_INs_slc32a1_or_gad_only)
granule.res.4.gruffi <- mature_INs_slc32a1_or_gad_only@misc$gruffi$'optimal.granule.res'
mature_INs_slc32a1_or_gad_only <- reassign.small.clusters(mature_INs_slc32a1_or_gad_only, ident = granule.res.4.gruffi) # will be stored in meta data column as "seurat_clusters.reassigned"
granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')

clustering = as.factor(granule.res.4.gruffi)

# ER stress 	GO:0034976
mature_INs_slc32a1_or_gad_only <- GO_score_evaluation(obj = mature_INs_slc32a1_or_gad_only, GO_term = go2, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
mature_INs_slc32a1_or_gad_only <- GO_score_evaluation(obj = mature_INs_slc32a1_or_gad_only, GO_term = go1, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Neurogenesis		GO:0022008
mature_INs_slc32a1_or_gad_only <- GO_score_evaluation(obj = mature_INs_slc32a1_or_gad_only, GO_term = go4, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)

# Create score names:
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i4 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go4))
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))

# Call Shiny app
mature_INs_slc32a1_or_gad_only <- Shiny.GO.thresh(obj = mature_INs_slc32a1_or_gad_only,
                                stress.ident1 = i2,
                                stress.ident2 = i1,
                                notstress.ident3 = i4,
                                plot.cluster.shiny = "seurat_clusters")

"Dont forget to click the button in the app: Save New Thresholds"

Seurat.utils::clUMAP('is.Stressed', label =T, obj = AD3)

cellIDs.keep <- which_names(!mature_INs_slc32a1_or_gad_only$'is.Stressed')
cellIDs.stressed <- which_names(mature_INs_slc32a1_or_gad_only$'is.Stressed')
mature_INs_stressed <- subset(x = mature_INs_slc32a1_or_gad_only, cells = cellIDs.stressed)  
mature_INs_mito_stress_filt_slc32a1_or_gad_only <- subset(x = mature_INs_slc32a1_or_gad_only, cells = cellIDs.keep)  

Seurat.utils::clUMAP('is.Stressed', label = F, obj = AD3_stressed)

#percentage of cells removed 
table(mature_INs_stressed@meta.data$orig.ident)

#total cells
table(mature_INs_mito_stress_filt@meta.data$orig.ident)

DimPlot(mature_INs_mito_stress_filt, reduction = "umap", group.by = "seurat_clusters")

save(mature_INs_mito_stress_filt_slc32a1_or_gad_only, file = "~/Desktop/mature_INs_mito_stress_filt_slc32a1_or_gad_only_20231003.RData")


# ######## PLOTS ########
# # stress assignment barplot cell numbers
# ggplot2::ggplot(combined.obj@meta.data, ggplot2::aes(x=combined.obj$orig.ident, fill=is.Stressed, stat = "count")) +
#   ggplot2::scale_fill_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) + ggplot2::theme_minimal() +
#   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   ggplot2::geom_bar() + ggplot2::xlab("") +
#   ggplot2::ggsave(paste0("Barplot_stressed.cells.png"), width = 3, height = 7)
#
# # stress assignment UMAP
# clUMAP(obj = combined.obj, ident = "is.Stressed", save.plot = F, label = F, legend = F)+#, splitby = "orig.ident") +
#   ggplot2::theme(text = ggplot2::element_text(size=20)) + Seurat::NoAxes() + ggplot2::ggtitle(ggplot2::element_blank())+ ggplot2::scale_color_manual(values=c(gg_color_hue(2)[2], gg_color_hue(2)[1])) +
#   ggplot2::ggsave(paste0("UMAP_Stress_cluster.png"), width = 10, height = 7)
