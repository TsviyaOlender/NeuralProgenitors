library(Seurat)
library(dplyr)

# Set working directory to the folder containing all 4 10x barcodes/features/matrices
setwd("Direcetory of the folder containing seurat objects")

# First, you can create Seurat objects from 10x files using the commented code below.
# The "folder" int he Read10x() function should be the directory of the folder containing barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
# raw.data <- Read10X(data.dir = "folder")
# SO <- CreateSeuratObject(counts = raw.data, project = "WT", assay = "RNA")

# read seurat objects
WT <- readRDS('WT')
KO <- readRDS('KO')

# Merge the Seurat objects into a single Seurat object
SO <- merge(x = WT, y = KO, add.cell.ids = c("WT", "KO"), project = "SO")
SO <- JoinLayers(SO)
# Change the data type of 'orig.ident' to factor for later processing
SO@meta.data$orig.ident <- as.factor(SO@meta.data$orig.ident)

# Filter out low-quality cells
SO[["percent.mt"]] <- PercentageFeatureSet(SO, pattern = "^MT-")
SO <- subset(SO, subset = nCount_RNA > 100 & nCount_RNA < 7.5e4 & percent.mt < 25 & percent.mt > 0.2)

# Normalization and cell cycle regression
all.genes <- rownames(SO)
SO <- NormalizeData(SO, normalization.method = "LogNormalize", scale.factor = 10000)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
SO <- CellCycleScoring(SO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SO <- ScaleData(SO, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)

# Save the cell cycle-regressed Seurat object into rds file
saveRDS(SO, "cell_cycle_regressed.rds")

