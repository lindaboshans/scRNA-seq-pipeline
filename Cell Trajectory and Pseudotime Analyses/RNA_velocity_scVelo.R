#RNA velocity with scVelo

install.packages("renv")
setwd("~/Desktop")
renv::init()
renv::install("reticulate")
library(reticulate)
reticulate::use_python("/usr/local/bin/python3.9", required = TRUE)
renv::use_python("/usr/local/bin/python3.9")


pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)


bioc_pkgs <- c(
  "SingleCellExperiment",
  "scater",
  "multtest"
)

#install in renv
renv::install(pkgs)

#install bioC packages
BiocManager::install(bioc_pkgs, update = FALSE)

#python packages
py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)

py_pkgs2 <- c(
  "anndata",
  "pandas",
  "velocyto"
)

#install in r
reticulate::py_install(py_pkgs2)

renv::snapshot()

# save metadata table:
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed$barcode <- colnames(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed)
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed$UMAP_1 <- mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed@reductions$umap@cell.embeddings[,1]
mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed$UMAP_2 <- mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed@reductions$umap@cell.embeddings[,2]
write.csv(gcdata_CC_regressed_CSS_integrated@meta.data, file='~/Desktop/gcdata_CC_regressed_CSS_integrated_metadata.csv', quote=F, row.names=F)
write.csv(counts, file='~/Desktop/gcdata_CC_regressed_CSS_integrated_RNA_counts.csv', quote=F, row.names=F)
write.table(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed@meta.data, file = "~/Desktop/mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_metadata.txt", quote=F, row.names=F, sep = '\t')

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0("~/Documents/", 'mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_counts.mtx'))
counts <- as.matrix(counts_matrix)

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed@reductions$umap@cell.embeddings, file='~/Desktop/mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_umap_coordinates.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='~/Desktop/mature_INs_mito_stress_filt_slc32a1_or_gad_only_CC_ribo_regressed_gene_names.csv',
  quote=F,row.names=F,col.names=F
)

####
####
####
# to use velocyto need to do the following in terminal;
#brew install llvm
#brew install boost
#brew install homebrew/science/hdf5 --enable-cxx
#mkdir -p ~/.R
#nano Makevars (to create file and paste in text below)
#CC=/usr/local/clang4/bin/clang
#CXX=/usr/local/clang4/bin/clang++
  #CXX11=/usr/local/clang4/bin/clang++
  #CXX14=/usr/local/clang4/bin/clang++
  #CXX17=/usr/local/clang4/bin/clang++
  #CXX1X=/usr/local/clang4/bin/clang++
  #LDFLAGS=-L/usr/local/clang4/lib
#then do below:
#Download clang-4.0.0 from http://r.research.att.com/libs/
#tar xvzf clang-4.0.0-darwin15.6-Release.tar in my Downloads directory
#cd usr/local/
  #sudo mv clang4 /usr/local/

library(Seurat)
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods", version = "3.12")
BiocManager::install("pcaMethods")
library(pcaMethods)
install.packages("magrittr")
devtools::install_github('velocyto-team/velocyto.R')
library(velocyto.R)
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
Sys.setenv('R_MAX_VSIZE'=24000000000)

###

run10x_48hr <- ReadVelocity(file = "~/Desktop/48hr.loom")
run10x_48hr_read_loom <- read.loom.matrices("~/Desktop/1148_GEX_48h1_48h2_hashed_3GEX_GRCh38-and-mm10.loom")
# Gather the spliced and unspliced estimates
emat <- run10x_48hr$spliced
nmat <- run10x_48hr$unspliced
emat <- emat[,colSums(emat)>=1e3]
# take embedding from the Seurat data object
# NOTE: This assumes you have a seurat data object loaded
# into R memory prior to using this script. STOP and rerun seurat
# pipeline if you do not have this loaded. In my case, my seurat object is simply myData
gc_48hr <- subset(gcdata, subset = timepoint == "48h_pD")
emb <- gc_48hr@reductions$umap@cell.embeddings
emb <- FetchData(object = gc_48hr, vars = c('UMAP_1','UMAP_2'))

# Estimate the cell-cell distances 
cell.dist <- as.dist(1-armaCor(t(emb)))
cell.dist <- as.dist(1-armaCor(t(FetchData(object = gc_48hr, vars = c('UMAP_1','UMAP_2')))))


#Filter genes based on the minimum average expresion magnitude (in at least one of the clusters), output total number of resulting valid genes
cluster.label <- gc_48hr@reductions$umap[[1]]
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
length(intersect(rownames(emat),rownames(emat)))

# This is a little bit of foo-magic that needs to be adjusted on a per-sample
# basis depending on the cell names and how you ran the pipeline. Each cell
# stored in the loom object and seurat have an ID, make sure these are the same.
# in this case -- i need to trim 28 characters off the front to match seurat object.
colnames(emat) <- paste(substring(colnames(emat),48,63),sep="")
colnames(nmat) <- paste(substring(colnames(nmat),48,63),sep="")
colnames(emat) <- paste0(colnames(emat),"_1")
colnames(nmat) <- paste0(colnames(nmat),"_1")
colnames(emat)
colnames(nmat)

# What this step does is essentially this:
# > head(colnames(emat))
# [1] "possorted_genome_bam_XL2S3:AAAGATGCATACTACGx"
# [2] "possorted_genome_bam_XL2S3:ACCTTTATCTTTAGTCx"
# [3] "possorted_genome_bam_XL2S3:AAGGAGCCACGCATCGx"

# > colnames(emat) <- paste(substring(colnames(emat),28,43),sep="")
# > head(colnames(emat))
# [1] "AAAGATGCATACTACG" "ACCTTTATCTTTAGTC" "AAGGAGCCACGCATCG" 

# Now the names in emat and nmat will match up to the cell names used in my seurat object

# I'm not sure what this parameter does to be honest. 0.02 default
# perform gamma fit on a top/bottom quantiles of expression magnitudes
fit.quantile <- 0.02

#need parallel forking to do below
install.packages("lme4",
                 repos=c("http://lme4.r-forge.r-project.org/repos",
                         getOption("repos")[["CRAN"]]))
library(lme4)

f <- function(i) {rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,
                                                              kCells=10,
                                                              cell.dist=cell.dist,
                                                              fit.quantile=fit.quantile,
                                                              n.cores=4)}

system.time(save1 <- lapply(1:100, f))

# Main velocity estimation
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,
                                            kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=4)


# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
gg <- TSNEPlot(myData)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)

p1 <- show.velocity.on.embedding.cor(emb,rvel.cd,n=30,scale='sqrt',
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=2,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                     do.par=F,cell.border.alpha = 0.1,
                                     n.cores=24,main="Cell Velocity")#,cc=p1$cc)


bm <- as.Seurat(x = run10x_48hr)
bm[["RNA"]] <- bm[["spliced"]]
#filter genes
counts <- GetAssayData(bm, slot="counts", assay="RNA")   
genes.percent.expression <- rowMeans(counts>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>1])  #select genes expressed in at least 1% of cells
counts.sub <- counts[genes.filter,]
filt_gc_48hr <- CreateSeuratObject(counts=counts.sub)
######

bm <- SCTransform(object = bm)
seurat_object <- Seurat::SCTransform(bm, assay = "spliced")
bm <- RunPCA(object = bm, verbose = FALSE)
bm <- FindNeighbors(object = bm, dims = 1:20)
bm <- FindClusters(object = bm)
bm <- RunUMAP(object = bm, dims = 1:20)
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
