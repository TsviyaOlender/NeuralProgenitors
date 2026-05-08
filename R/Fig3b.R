library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)

# load cell cycle regressed Seurat object and subset rostral samples
setwd("Direcetory of the folder containing seurat objects")
SO <- readRDS("cell_cycle_regressed.rds")
SO <- SetIdent(SO, value = 'orig.ident')

# Perform UMAP and clustering
SO <- FindVariableFeatures(SO)
SO <- RunPCA(SO, features = VariableFeatures(object = SO))
SO <- RunUMAP(SO, dims = 1:10, reduction = "pca", return.model = TRUE)
SO <- FindNeighbors(SO, dims = 1:10)
SO <- FindClusters(SO, resolution = 1.1) # Resolution of 1.1 is for sample data. Use resolution of 0.1 for the full data. 

# Fig 3b
SO <- SetIdent(SO, value = 'orig.ident')
DimPlot_scCustom(subset(SO, ident = c('WT')), group.by = 'seurat_clusters', reduction = 'umap', label = FALSE, colors_use = c('#406B80', "#BE1E2D", '#102B40', '#A3C2CC'), pt.size =  4) + NoLegend() +NoAxes()+ ggtitle(NULL) # 2000x1600
DimPlot_scCustom(subset(SO, ident = c('KO')), group.by = 'seurat_clusters', reduction = 'umap', label = FALSE, colors_use = c('#406B80', "#BE1E2D", '#102B40', '#A3C2CC'), pt.size =  4) + NoLegend() +NoAxes()+ ggtitle(NULL) # 2000x1600

# Feature plots in Supplementary Fig S3
marker.genes = c('NRP2', 'PRKCB', 'SNAI1', 'DDIT3', 'FGF8', 'VEGFA', 'STAT3', 'PDGFA', 'ACVR1', 'PPP1R14A', 'ARHGEF4', 'NFKBIB')
marker.genes %in% row.names(SO)

saving_dir <- 'directory to save'
for (gene in marker.genes){
  directory <- paste(saving_dir, gene, '.png', sep = '')
  print(directory)
  FeaturePlot_scCustom(rostral, features = c(gene), order = TRUE, pt.size = 0.4, slot = 'data', na_cutoff = 0.01, split.by = 'orig.ident')  + NoAxes()
  ggsave(directory, width = 2000, height = 1000, units = "px")
}

saveRDS(SO, 'clustered.rds')

