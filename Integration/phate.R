install.packages("phateR")
library(reticulate)
reticulate::py_install("phate", pip=TRUE)
py_config()
reticulate::py_discover_config(required_module = "phate")
reticulate::import("phate")
library(phateR)

### handles the Seurat to phate conversio

### first, grab the input required for phate (here we are using the normalized data stored in Seurat
seurat_data <- as.data.frame(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed@assays$RNA@data)

## reshape for input into PHATE
phate_data_input <- t(seurat_data)

## run phate -- this is a default run, feel free to tune params
phate_output <- phate(phate_data_input, ndim=30) 
phate_output <- phate(phate_data_input, t=10, ndim=30, init=phate_output) 

## quick sanity check of the phate embeddings
ggplot(phate_output, aes(x=PHATE1, y=PHATE2)) +
  geom_point()

## stash the embeddings back into the seurat object as a dimension reduction object
gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed[["PHATE"]] <- CreateDimReducObject(embeddings = phate_output$embedding, key = "PHATE_", assay = DefaultAssay(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed))

## plot using seurat's default tools
DimPlot(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed , reduction = "PHATE", group.by = "timepoint")

save(gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed, file = "~/Desktop/gcdata_mito_filtered_all_10_wk2_15_stress_0.56_CC_regressed_w_default_phate.RData")

