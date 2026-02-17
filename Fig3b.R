library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)

# load cell cycle regressed Seurat object and subset rostral samples
setwd("your working directory")
SO <- readRDS("cell_cycle_regressed.rds")
SO <- SetIdent(SO, value = 'orig.ident')
rostral <- subset(x = NDE1, ident = c('KOR', 'WTR'))
rm(NDE1)

# Perform UMAP and clustering
rostral <- RunPCA(rostral, features = VariableFeatures(object = rostral))
rostral <- RunUMAP(rostral, dims = 1:10, reduction = "pca", return.model = TRUE)
rostral <- FindNeighbors(rostral, dims = 1:10)
rostral <- FindClusters(rostral, resolution = 0.1)

# Fig 3b
rostral <- SetIdent(rostral, value = 'orig.ident')
DimPlot_scCustom(subset(rostral, ident = c('WTR')), group.by = 'seurat_clusters', reduction = 'umap', label = FALSE, colors_use = c('#406B80', "#BE1E2D", '#102B40', '#A3C2CC'), pt.size =  4) + NoLegend() +NoAxes()+ ggtitle(NULL) # 2000x1600
DimPlot_scCustom(subset(rostral, ident = c('KOR')), group.by = 'seurat_clusters', reduction = 'umap', label = FALSE, colors_use = c('#406B80', "#BE1E2D", '#102B40', '#A3C2CC'), pt.size =  4) + NoLegend() +NoAxes()+ ggtitle(NULL) # 2000x1600

# Feature plots in Supplementary Fig S3
marker.genes = c('NRP2', 'PRKCB', 'SNAI1', 'DDIT3', 'FGF8', 'VEGFA', 'STAT3', 'PDGFA', 'ACVR1', 'PPP1R14A', 'ARHGEF4', 'NFKBIB')
marker.genes %in% row.names(rostral)

rostral <- SetIdent(rostral , value = 'orig.ident')
WTR <- subset(rostral, ident = c('WTR'))
KOR <- subset(rostral, ident = c('KOR'))

orident <- rostral$orig.ident
orident[which(rostral$orig.ident == 'WTR')] <- 'WT'
orident[which(rostral$orig.ident == 'KOR')] <- 'KO'
orident <- factor(orident, levels = c('WT', 'KO'))
rostral$orig.ident <- orident


saving_dir <- 'directory to save'
for (gene in marker.genes){
  directory <- paste(saving_dir, gene, '.png', sep = '')
  print(directory)
  FeaturePlot_scCustom(rostral, features = c(gene), order = TRUE, pt.size = 0.4, slot = 'data', na_cutoff = 0.01, split.by = 'orig.ident')  + NoAxes()
  ggsave(directory, width = 2000, height = 1000, units = "px")
}

saveRDS(rostral, 'rostral.rds')

