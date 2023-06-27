# Integration of two groups of PBMC (control and stimulated cells)
# The response to interferon caused cell type specific gene expression changes that makes a joint analysis of all the data difficult
# see https://satijalab.org/seurat/archive/v3.2/immune_alignment.html


library(Seurat)
library (SeuratData)

# InstallData ("ifnb")
data("ifnb")
ifnb <- UpdateSeuratObject(ifnb)
print(ifnb)


ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list
# CTRL and STIM list


ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Anchors between the two datasets
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
# Merging/integrating dataset 1 into 2
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


# run a single integrated analysis on all types of cells
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
library (cowplot)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")


# Make a column summarizing the celltype (or Seurat cluster) and the treatment  
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
table (immune.combined@meta.data$celltype.stim)
table (immune.combined@meta.data$seurat_annotations)


# Differential expression in a particular cell type, here CD8 T cells
mycl <- immune.combined@meta.data
mycl <- mycl[mycl$seurat_annotations == "CD8 T", ]
table (mycl$seurat_clusters)
# equivalent to cluster 5

## deprecated for Presto. Do not use
#b.interferon.response <- FindMarkers(immune.combined, ident.1 = "5_STIM", ident.2 = "5_CTRL", verbose = FALSE)
#head(b.interferon.response, n = 15)



#### Presto performs a fast Wilcoxon rank sum test and auROC analysis
# See https://github.com/immunogenomics/presto/blob/master/vignettes/getting-started.Rmd

#install_github('immunogenomics/presto')

library(presto)

# here, we a comparing 5_STIM and 5_CTRL groups, however, both directions are reported !
res_AB <- wilcoxauc (immune.combined, "celltype.stim", assay = 'data', groups_use = c('5_STIM', '5_CTRL'))
res_AB <- res_AB[order (res_AB$padj), ]
res_AB <- res_AB[res_AB$group == "5_STIM", ]
head (res_AB)



#### fgsea

library(msigdbr)
library(fgsea)
library(dplyr)
library (tidyverse)
library(ggplot2)

msigdbr_species()

m_df <- msigdbr(species = "Homo sapiens", category = "C5")
head(m_df)

fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
head (fgsea_sets)


## select only the feature and auc columns for fgsea, which statistics to use is an open question !

genes <- res_AB %>%
  				  dplyr::filter(group == "5_STIM") %>%
                  arrange(desc(auc)) %>% 
                  dplyr::select(feature, auc)

head (genes)

ranks <- deframe(genes)
head(ranks)


fgseaRes <- fgseaMultilevel (fgsea_sets, stats = ranks, minSize=15, maxSize=500)

fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
head (data.frame (fgseaResTidy))


plotEnrichment(fgsea_sets[["GOBP_RESPONSE_TO_VIRUS"]], ranks)

# The Enrichement Score (ES) is calcuated by some metric that ES is positive if the gene set is located in the top of the pre-ranked gene list. 
# The ES is negative if the gene set is located in the bottom of the pre-ranked gene list.



