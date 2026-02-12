library(SingleR)
library(SingleCellExperiment)
orly<-readRDS("/home/labs/olenderlab/lvzvia/reinerlab/Tamar_scRNA_NDE_jan23/analysis/seurat.analysis_30Jan23_25PCA_res0.6.rds")
orly<-NormalizeData(orly)
orly<-FindVariableFeatures(orly)
orly<-ScaleData(orly)
DefaultAssay(orly)<-"RNA"
#
sub<-readRDS("/home/labs/olenderlab/lvzvia/MyRScripts/scRNA_pipline/human_embryonic_atlas/ClassRadialGlia_v3_4Aug25_sketch_100K.rds")
DefaultAssay(sub)<-'sketch'
counts_mat <- GetAssayData(sub, assay = "sketch", slot = "counts")
query <- CreateSeuratObject(counts = counts_mat, meta.data = sub@meta.data)
#
query<-NormalizeData(query)
query <- FindVariableFeatures(query, verbose = FALSE)
query <- ScaleData(query)
query <- RunPCA(query)
query<-RunUMAP(query, dims = 1:40, verbose = FALSE)
#
#flex_xen_common_genes <- intersect(rownames(query), rownames(orly))

pancreas.anchors <- FindTransferAnchors(reference = query, query = orly, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = query$Region, dims = 1:30)
orly <- AddMetaData(orly, metadata = predictions)
Idents(orly)<-orly$predicted.id
DimPlot_scCustom(orly, split.by = "genotypes",label = F)
FeaturePlot_scCustom(orly,features = "prediction.score.max",split.by = "genotypes")

Cluster_Highlight_Plot(seurat_object = orly, cluster_name = "forebrain", highlight_color = "forestgreen",
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
orly$genotypes <- factor(orly$genotypes, levels = desired_order)
DimPlot_scCustom(orly,colors_use = my_colors,split.by = "genotypes",label=F,pt.size=1)

svg("UMAP_by_Miri_1.svg", width = 15, height = 8)  # size in inches
DimPlot_scCustom(orly,colors_use = my_colors,split.by = "genotypes",label=F,pt.size=2)
dev.off()
