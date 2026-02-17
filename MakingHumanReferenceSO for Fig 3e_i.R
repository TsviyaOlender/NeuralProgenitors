library(Seurat)
library(hdf5r)
library(BiocManager)
BiocManager::install("rhdf5")
library(rhdf5)

# Reading in human embryonic brain atlas
# The original work was:Braun, Emelie, et al. "Comprehensive cell atlas of the first-trimester developing human brain." Science 382.6667 (2023): eadf1226.
# The data was retrieved from: https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077
setwd("your working directory")
h5ls('HumanFetalBrainPool.h5')

# Retrieve each fields
age <- h5read(file = h5_file_path, name = "shoji/Age")
AnnotationName <- h5read(file = h5_file_path, name = "shoji/AnnotationName")
CellClass <- h5read(file = h5_file_path, name = "shoji/CellClass")
Class <- h5read(file = h5_file_path, name = "shoji/Class")
Region	 <- h5read(file = h5_file_path, name = "shoji/Region")
Tissue	 <- h5read(file = h5_file_path, name = "shoji/Tissue")
CellID	<- h5read(file = h5_file_path, name = "shoji/CellID")
Gene	<- h5read(file = h5_file_path, name = "shoji/Gene")
Subdivision	<- h5read(file = h5_file_path, name = "shoji/Subdivision")
Subregion	<- h5read(file = h5_file_path, name = "shoji/Subregion")
ValidCells	<- h5read(file = h5_file_path, name = "shoji/ValidCells")
ValidGenes	<- h5read(file = h5_file_path, name = "shoji/ValidGenes")

# generate seurat object
valid.gene.idx <- which(ValidGenes == TRUE)
idx <- which( (age == 7) & CellClass == 'Radial glia')
Expression	<- h5read(file = h5_file_path, name = "shoji/Expression", index = list(valid.gene.idx, idx)) # genes cells
rownames(Expression) <- Gene[valid.gene.idx]
colnames(Expression) <- CellID[idx]
Expression <- Expression[which(!duplicated(rownames(Expression))),]
meta.data <- data.frame(age[idx], CellClass[idx], Region[idx], Tissue[idx], Subdivision[idx], Subregion[idx])
rownames(meta.data) <- CellID[idx]
colnames(meta.data) <- c('age', 'CellClass', 'Region','Tissue', 'Subdivision', 'Subregion')
SO <- CreateSeuratObject(counts = Expression, meta.data = meta.data)
saveRDS(SO, 'HumanReference.rds')
