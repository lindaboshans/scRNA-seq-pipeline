
library(RRHO2)


wk4_map2_slc32a1_only <- AddMetaData(wk4_map2_slc32a1_only, metadata = "wk4", col.name = "structure")

wk4_day2 <- merge(wk4_map2_slc32a1_only, y = day2, add.cell.ids = c("Wk4", "day2"), project = "DE")
merge_stri_removed <- merge(wk4_map2_slc32a1_only, y = GW20_brain_regions_combined_cells_filtered_neurons_GADexp_intneurons_stri_removed, add.cell.ids = c("Wk4", "GW20_no_stri"), project = "mature_INs")


GW20_brain_regions_combined_cells_filtered_neurons_GADexp_intneurons_hypo_removed <- subset(GW20_brain_regions_combined_cells_filtered_neurons_GADexp_intneurons, subset = structure == "hypothalamus", invert = TRUE)
GW20_brain_regions_combined_cells_filtered_neurons_GADexp_intneurons_stri_removed <- subset(GW20_brain_regions_combined_cells_filtered_neurons_GADexp_intneurons, subset = structure == "striatum", invert = TRUE)

Idents(merge_hypo_removed) <- "structure"
Idents(GW20_brain_regions_combined_cells_filtered) <- "structure"
Idents(wk4_day2) <- "timepoint"

g <- FindMarkers(wk4_day2, ident.1 = '4wk_glia')
write.csv(g, file = "~/Documents/Ruiqi_scRNAseq/wk4_vs_day2_findallmarkers_DE_genes_for_RRHO2_20231106.csv")

g <- FindAllMarkers(GW20_brain_regions_combined_cells_filtered_neurons)
write.csv(g, file = "~/Documents/Ruiqi_scRNAseq/wk4_vs_gw20_findallmarkers_DE_genes_stri_removed_for_RRHO2_20231103.csv")

f <- FindAllMarkers(GW20_brain_regions_combined_cells_filtered)
write.csv(f, file = "~/Documents/Ruiqi_scRNAseq/GW20_all_cells_findallmarkers_DE_genes_for_RRHO2_20231103.csv")


class(GW20_brain_regions_combined_cells_filtered_neurons@assays$RNA@counts)


set.seed(15213)
nGenes <- 2000
nDE <- 200
Genes <- paste0("Genes",1:nGenes)

## For up-regulated genes, the input score should be calculated using-log10(pvalue) * 1;
## For down-regulated genes, the input score should be calculated using-log10(pvalue) * (-1);

list1_pvalue_1_200 <- runif(nDE,0,0.05) ## up-regulated genes
list1_pvalue_201_400 <- runif(nDE,0,0.05) ## down-regulated genes
list1_pvalue_401_2000 <- runif(nGenes - 2 * nDE,0,1) ## non-changed genes
list1_DDE <- c(-log10(list1_pvalue_1_200), 
               -log10(list1_pvalue_201_400) * (-1), 
               -log10(list1_pvalue_401_2000) * sample(c(1,-1), length(list1_pvalue_401_2000), replace = TRUE))

gene_list1 <- data.frame(Genes=Genes,DDE = list1_DDE, stringsAsFactors = FALSE)
gene_list1 <- data.frame(Genes=gene_list1$gene,DDE = gene_list1$DDE, stringsAsFactors = FALSE)
gene_list1 <- read.csv("~/Documents/Ruiqi_scRNAseq/wk4 map2 slc32a1 RRHO2 input DE vs day2.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)


list2_pvalue_1_200 <- runif(nDE,0,0.05)
list2_pvalue_201_400 <- runif(nDE,0,0.05) 
list2_pvalue_401_2000 <- runif(nGenes - 2 * nDE,0,1)
list2_DDE <- c(-log10(list2_pvalue_1_200), -log10(list2_pvalue_201_400) * (-1), 
               -log10(list2_pvalue_401_2000) * sample(c(1,-1), length(list2_pvalue_401_2000), replace = TRUE))

gene_list2 <- data.frame(Genes=Genes,DDE = list2_DDE, stringsAsFactors = FALSE)
gene_list2 <- data.frame(Genes=gene_list2$gene,DDE = gene_list2$DDE, stringsAsFactors = FALSE)
gene_list2 <- read.csv("~/Documents/Ruiqi_scRNAseq/hypo RRH02 input final GW20 filtered neurons.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)
head(gene_list2)

RRHO_obj <-  RRHO2_initialize(gene_list1, gene_list2, labels = c("gene_list1", "gene_list2"), log10.ind=TRUE)

RRHO2_heatmap(RRHO_obj)
