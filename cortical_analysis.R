library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(rio)
library(grid)
library(gridExtra)
library(scCustomize)
##############################################################################
# FUNCTION: filter seurat object
###############################################################################
filter_seurat_object_1<- function(obj.seurat,mitho_max,genes_min,genes_max,umi_min,umi_max,obj.num) {
  obj.seurat[["percent.mt"]] <- PercentageFeatureSet(obj.seurat , pattern = "^MT-")
  if(obj.num > 1){
    p.before<-VlnPlot(obj.seurat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),split.by="orig.ident", ncol = 3,pt.size = 0.05)
    print(p.before)
  }else{
    p.before<-VlnPlot(obj.seurat, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.05)
    print(p.before)
  }
  obj.seurat.f <- subset(obj.seurat, subset = percent.mt < mitho_max &
                           nFeature_RNA>genes_min & nFeature_RNA < genes_max &
                           nCount_RNA > umi_min & nCount_RNA < umi_max)
  # plot after filtration
  obj.seurat.f[["percent.mt"]] <- PercentageFeatureSet(obj.seurat.f , pattern = "^MT-")
  if(obj.num > 1){
    p.after<-VlnPlot(obj.seurat.f, features = c("nFeature_RNA", "nCount_RNA","percent.mt"),split.by="orig.ident", ncol = 3,pt.size = 0.05)
    print(p.after)
  }else{
    p.after<-VlnPlot(obj.seurat.f, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.05)
    print(p.after)
  }
  ########
  # print parameters
  filter.parameters <- matrix(c(rep(0,5)), ncol=5, byrow=TRUE)
  filter.parameters[1,]=c(mitho_max,genes_min,genes_max,umi_min,umi_max)
  colnames(filter.parameters) <- c("mitho_max","genes_min","genes_max","umi_min","umi_max")
  print(filter.parameters)
  
  return(obj.seurat.f)
}
##############################################################################
# FUNCTIONS end
#############################################################################
# input parameters
data.folder="/home/labs/olenderlab/lvzvia/reinerlab/Tamar_scRNA_NDE_jan23/cell_ranger_analysis_NDE_KO/outs/per_sample_outs"
save.dir = "/home/labs/olenderlab/lvzvia/reinerlab/Tamar_scRNA_NDE_jan23/analysis"
experiments = c("NDE1_KO_3","NDE1_KO_4","W3_E3","W3_EP3")
experimentNames = c("NDE1_KO_3","NDE1_KO_4","W3_E3","W3_EP3")
samples.num=4
pca.num=25 # Num of PCA
res.num = 0.6 # resolution
# seurat out file
seurat_out_file="seurat.analysis_corr_Jan23.rds"
############################################################################################
# load data and create seurat objects
############################################################################################
for (i in 1:samples.num){
  experimentLocation = paste(data.folder,"/",experiments[i],"/count/sample_filtered_feature_bc_matrix",sep="")
  experiment.data <- Read10X(data.dir =experimentLocation);
 seurat.obj<- CreateSeuratObject(counts  = experiment.data$`Gene Expression`, 
                                 project = experimentNames[i], min.cells = 3, min.features = 0)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj , pattern = "^MT-")
  assign(experimentNames[i],seurat.obj)
  
}

##############################################################################
# filter each individual sample
#NDE1_KO_3","NDE1_KO_4","W3_E3","W3_EP3"
#############################################################################
NDE1_KO_3.f<-filter_seurat_object_1(NDE1_KO_3,10,100,10000,200,30000,1)
NDE1_KO_4.f<-filter_seurat_object_1(NDE1_KO_4,10,100,10000,200,30000,1)
W3_E3.f<-filter_seurat_object_1(W3_E3,10,0,10000,200,30000,1)
W3_EP3.f<-filter_seurat_object_1(W3_EP3,10,0,10000,200,30000,1)
seurat.obj <- merge(x=NDE1_KO_3,y=list(NDE1_KO_4,W3_E3,W3_EP3),
                      add.cell.ids = experimentNames)
table(seurat.obj$orig.ident)
VlnPlot(seurat.obj  , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),split.by="orig.ident", ncol = 3,pt.size = 0.05)
seurat.obj.f <- merge(x=NDE1_KO_3.f,y=list(NDE1_KO_4.f,W3_E3.f,W3_EP3.f),
                      add.cell.ids = experimentNames)
VlnPlot(seurat.obj.f  , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),split.by="orig.ident", ncol = 3,pt.size = 0.05)
table(seurat.obj.f$orig.ident)
####################################################################################
# rename
#NDE1_KO_3 NDE1_KO_4     W3_E3    W3_EP3 
#1250      3308      1026      2305 
#NDE1_KO_4 is a WT, and that 'W3_E3' is KO. 

seurat.obj.f$corrected.ident<-c(rep("KO_3",1250),rep("W3_e3",3308),rep("KO_4",1026),rep("W3_ep3",2305))
seurat.obj.f$genotypes<-c(rep("KO",1250),rep("WT",3308),rep("KO",1026),rep("WT",2305))
#######################################################################################
# RNA normalization, calculate cell cycle score
#######################################################################################
seurat.obj.f <- NormalizeData(seurat.obj.f)
seurat.obj.f <- FindVariableFeatures(seurat.obj.f, selection.method = "vst")
seurat.obj.f <- ScaleData(seurat.obj.f, features = row.names(seurat.obj.f))
#
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat.obj.f <- CellCycleScoring(seurat.obj.f, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#
#############################################################################
# SCT Normalization
####################################################################################
seurat.obj.f <- SCTransform(seurat.obj.f, verbose = FALSE,vst.flavor = "v2")
# PCA, UMAP with object name
seurat.obj.f <- RunPCA(seurat.obj.f, verbose = FALSE)
ElbowPlot(seurat.obj.f,ndims=40)
print(seurat.obj.f[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat.obj.f, dims = 1:2, reduction = "pca")
DimHeatmap(seurat.obj.f, dims = 20:30, cells = 500, balanced = TRUE)
FeaturePlot(seurat.obj.f,features=c("NDE1"),reduction="umap",split.by = "orig.ident")
############################################################################
# clustering and dimentionality reduction
seurat.obj.f <- FindNeighbors(seurat.obj.f, dims = 1:pca.num)
seurat.obj.f <- FindClusters(seurat.obj.f, resolution = res.num)
seurat.obj.f <- RunUMAP(seurat.obj.f, reduction = "pca", dims = 1:pca.num, verbose = FALSE)
p1<-DimPlot(seurat.obj.f, reduction = "umap",label=TRUE,label.size = 4,split.by="genotypes")
p2<-DimPlot(seurat.obj.f, reduction = "umap",label=TRUE,label.size = 6)
p1+p2

UMAP.cluster.file.bySample<-paste(save.dir,"/","UMAP_clusters_",pca.num,"PCA_res_",res.num,"_bySample.png",sep="")
UMAP.clusters.file = paste(save.dir,"/","UMAP_clusters_",pca.num,"PCA_res_",res.num,".png",sep="")
ggplot2::ggsave(filename = UMAP.cluster.file.bySample, plot = p1) 
ggplot2::ggsave(filename = UMAP.clusters.file, plot = p2) 

############################################################################
# Get number of cells per cluster and per sample of origin
#############################################################################
cell.count<-table(seurat.obj.f@meta.data[["SCT_snn_res.0.6"]], seurat.obj.f@meta.data$genotypes)
# Get number of cells per cluster and per sample of origin

clusters<-as.integer(levels(seurat.obj.f@meta.data[["SCT_snn_res.0.6"]]))
cluster.num=max(clusters)+1
cell.count.norm <- matrix(c(rep(0,cluster.num*samples.num)), ncol=samples.num, byrow=TRUE)

colnames(cell.count.norm) <- c("%KO","%WT")
rownames(cell.count.norm) <- c(0:(cluster.num-1))
cell.count.norm <- as.table(cell.count.norm)
for (i in 1:samples.num){
  cell.count.norm[,i]<-(cell.count[,i]/sum(cell.count[,i]))*100
}

final_count<-cbind(cell.count,cell.count.norm)
final.count.f<-as.data.frame(final_count)
final.count.f$cluster<-c(0:max(clusters))
cell.count.file = paste(save.dir,"/","cell_count_",pca.num,"PCA_res_",res.num,".xlsx",sep="")
rio::export(final.count.f,file=cell.count.file,overwrite=TRUE )
# plot
cell.count.norm<-t(cell.count.norm)
barplot(cell.count.norm,col = c("#1b98e0", "#353436"),beside = TRUE,
        xlab="cluster",ylab="% of cells",main="%cells per cluster",legend.text=TRUE)
###############################################################
# differential expression analysis
###############################################################
seurat.obj.f.markers <- FindAllMarkers(seurat.obj.f, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.5)
#TOP10
seurat.obj.f.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#top2
seurat.obj.f.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top2
HP2<-DoHeatmap(seurat.obj.f, features = c(top2$gene)) + NoLegend()
# save. file names
markers.file = paste(save.dir,"/","markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
top2.file = paste(save.dir,"/","top2..markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
top10.file = paste(save.dir,"/","top10.markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
top2.HM.file = paste(save.dir,"/","top2.heatmap_",pca.num,"PCA_res_",res.num,".png",sep="")
# write
write.csv(seurat.obj.f.markers,file=markers.file)
write.csv(top2,file=top2.file)
write.csv(top10,file=top10.file)
ggplot2::ggsave(filename = top2.HM.file, plot = HP2) 
############################################################################
# Test DE control versus KO
# across clusters
#############################################################
seurat.obj.f$cluster.genotype <- paste(seurat.obj.f$seurat_clusters, seurat.obj.f$genotypes, sep="_")
Idents (seurat.obj.f) <- seurat.obj.f$cluster.genotype
report<-""

#need to see which group has fewer than 3 cells
table(Idents (seurat.obj.f))
#clusters 0,1,2 AND 6 has less than 3 cells in one of the groups
clusters.to.check<-c(0,5,6,7,8,9,10,11,12)
for (j in c(1:9))
{  #number of clusters
  ident1 <- paste0(clusters.to.check[j],"_KO")
  ident2 <- paste0(clusters.to.check[j],"_WT")
   diffgenes.cluster <- FindMarkers(seurat.obj.f, ident.1 = c(ident1), ident.2=c(ident2),
                                   logfc.threshold=0.2, min.pct = 0.2 ,only.pos = FALSE)
  
  if (nrow(diffgenes.cluster)>0){ #there are diff genes
    diffgenes.cluster$contrast<-paste0(ident1, "_vs_", ident2)
    diffgenes.cluster$gene<- rownames(diffgenes.cluster)
    print (paste0("Analysis of cluster ", j))
    head(diffgenes.cluster)
    report<-rbind(report,diffgenes.cluster)
  }
}
write.table(report , file = "Per_cluster_Ctrl_vs_Mutant_30PCA_res0.6.def.nointegration.txt", sep = "\t",row.names = TRUE)
Idents (seurat.obj.f) <- seurat.obj.f$seurat_clusters
#########################################################################
# save data
##########################################################################
object.file.name= paste(save.dir,"/",seurat_out_file,sep="")
saveRDS(seurat.obj.f, file = object.file.name)
sessionInfo() %>% capture.output(file="session_info.txt")
############################################################