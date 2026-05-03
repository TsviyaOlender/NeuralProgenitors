library(SingleR)
library(SingleCellExperiment)
my.obj<-readRDS("seurat_obj_100cells.rds")
my.obj<-NormalizeData(my.obj)
my.obj<-FindVariableFeatures(my.obj)
my.obj<-ScaleData(my.obj)
DefaultAssay(my.obj)<-"RNA"
#
sub<-readRDS("/home/labs/olenderlab/lvzvia/MyRScripts/scRNA_pipline/human_embryonic_atlas/ClassRadialGlia_v3_4Aug25.rds")
counts_mat <- GetAssayData(sub, assay = "RNA", slot = "counts")
query <- CreateSeuratObject(counts = counts_mat, meta.data = sub@meta.data)
#
query<-NormalizeData(query)
query <- FindVariableFeatures(query, verbose = FALSE)
query <- ScaleData(query)
query <- RunPCA(query)
query<-RunUMAP(query, dims = 1:40, verbose = FALSE)
#
proj.anchors <- FindTransferAnchors(reference = query, query = my.obj, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = proj.anchors, refdata = query$Region, dims = 1:30)
my.obj <- AddMetaData(my.obj, metadata = predictions)
Idents(my.obj)<-my.obj$predicted.id
DimPlot_scCustom(my.obj, split.by = "genotypes",label = F)
FeaturePlot_scCustom(my.obj,features = "prediction.score.max",split.by = "genotypes")

Cluster_Highlight_Plot(seurat_object = my.obj, cluster_name = "forebrain", highlight_color = "forestgreen",
                       background_color = "lightgray")
my_colors=c("#5A5156FF", "#E4E1E3FF", "#F6222EFF" ,"#FE00FAFF" ,"#16FF32FF", "#3283FEFF", "#FEAF16FF", "#B00068FF")
my_colors <- c(
  "#FF4E4E",  # Bright Red
  "#4EA5FF",  # Sky Blue
  "#63E05B",  # Lime Green
  "#C17DFF",  # Lavender
  "#FFA742",  # Bright Orange
  "#5A5156FF",  # Light Yellow
  "#D98A5E",  # Tan/Bright Brown
  "#FFB3E6"   # Light Pink
)
desired_order <- c("WT","KO")
my.obj$genotypes <- factor(my.obj$genotypes, levels = desired_order)
DimPlot_scCustom(my.obj,colors_use = my_colors,split.by = "genotypes",label=F,pt.size=1)

svg("UMAP_by_Miri_1.svg", width = 15, height = 8)  # size in inches
DimPlot_scCustom(my.obj,colors_use = my_colors,split.by = "genotypes",label=F,pt.size=2)
dev.off()
