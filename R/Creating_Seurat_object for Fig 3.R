library(Seurat)
library(dplyr)

# Set working directory to the folder containing all 4 10x barcodes/features/matrices
setwd("your working directory")

# Read in each 10x files and create Seurat objects
WTR_raw <- Read10X(data.dir = "WTR")
WTR <- CreateSeuratObject(counts = WTR_raw, project = "WTR", assay = "RNA")

WTC_raw <- Read10X(data.dir = "WTC")
WTC <- CreateSeuratObject(counts = WTC_raw, project = "WTC", assay = "RNA")

KOR_raw <- Read10X(data.dir = "KOR")
KOR <- CreateSeuratObject(counts = KOR_raw, project = "KOR", assay = "RNA")

KOC_raw <- Read10X(data.dir = "KOC")
KOC <- CreateSeuratObject(counts = KOC_raw, project = "KOC", assay = "RNA")

# Merge the Seurat objects into a single Seurat object
SO <- merge(x = WTR, y = c(WTC, KOR, KOC), add.cell.ids = c("WTR", "WTC", "KOR", "KOC"), project = "SO")

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

