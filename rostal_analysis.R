library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(rio)
library(grid)
library(gridExtra)
library(scCustomize)
####################################################################################
# FUNCTION::: Remove genes that affect the PCA and invovled in Cell Cycle
###################################################################################
remove_genes_from_matrix<-function(the.object,genes.file){
  cc.genes<-readLines(genes.file)
  all.features = rownames(x = the.object)
  rowNum<-c()
  gene.nums<-as.integer(length(cc.genes))
  for(i in 1:gene.nums){
    rowNum[i]<-which(all.features ==cc.genes[i])
  
  }
  all.features <- all.features[-rowNum]

  the.object<-subset(x = the.object, features =all.features )
  return(the.object)
}
####################################################################################
# FUNCTION end
###################################################################################
# input parameters
###################################################################################
#cellranger files
experiments = c("Sample_8652-JB-1-GEX_GCAGTATAGG-GTGCACGGAA",
                "Sample_8652-JB-3-GEX_CCCACCACAA-AAGCGGAGGT")
experimentNames = c("WR","KOR")
data.folder="cellranger"
save.dir = "/home/labs/olenderlab/lvzvia/reinerlab/neuromorpho/Aug23"
samples.num=2
pca.num=30 # Num of PCA
res.num = 0.4 # resolution
##
# seurat out file
seurat_out_file="cleanSeuratobj_3Aug23.rds"
############################################################################################
# load data and create seurat objects
for (i in 1:samples.num){
  experimentLocation = paste(data.folder,"/",experiments[i],"/filtered_feature_bc_matrix",sep="")
  experiment.data <- Read10X(data.dir =experimentLocation);
  seurat.obj<- CreateSeuratObject(counts  = experiment.data, project = experimentNames[i], min.cells = 3, min.features = 0)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj , pattern = "^MT-")
  assign(experimentNames[i],seurat.obj)
  
}
seurat.obj<- merge(x=WR,y=list(KOR),add.cell.ids = experimentNames)
VlnPlot(seurat.obj,features = c("nFeature_RNA","nCount_RNA","percent.mt"))
seurat.obj.f <- subset(seurat.obj, subset = percent.mt <10 &
                         nCount_RNA > 200 & nCount_RNA < 50000&
                       nFeature_RNA > 1000 & nFeature_RNA < 10000)
VlnPlot(seurat.obj.f,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0)
#######################################################################################
# Normalize and calculate cell cycle score
#######################################################################################
seurat.obj.f <- NormalizeData(seurat.obj.f)
seurat.obj.f <- FindVariableFeatures(seurat.obj.f, selection.method = "vst")
seurat.obj.f <- ScaleData(seurat.obj.f, features = row.names(seurat.obj.f))
#
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat.obj.f <- CellCycleScoring(seurat.obj.f, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#
seurat.obj.f <- SCTransform(seurat.obj.f, vst.flavor = "v2", vars.to.regress = c("S.Score", "G2M.Score"))
seurat.obj.f <- RunPCA(seurat.obj.f, features = c(s.genes, g2m.genes),approx=FALSE)
DimPlot(seurat.obj.f)
####################################################################################
# Remove genes that affect the PCA and invovled in Cell Cycle
###################################################################################
new.object<-remove_genes_from_matrix(new.object,"cell_division_genes.F.txt")
new.object <- ScaleData(new.object, features=row.names(new.object))
new.object <- RunPCA(new.object, approx=FALSE)
new.object<-RunUMAP(new.object,reduction = "pca", dims = 1:pca.num, verbose = FALSE)
DimPlot_scCustom(new.object,split.by = 'orig.ident')
################################################################################
# clustering
####################################################################################
seurat.obj.f <- FindNeighbors(seurat.obj.f, dims = 1:pca.num)
seurat.obj.f <- FindClusters(seurat.obj.f, resolution = res.num)
DimPlot_scCustom(seurat.obj.f,split.by = 'orig.ident')
###############################################################
# differential expression analysis
###############################################################
seurat.obj.markers <- FindAllMarkers(seurat.obj.f, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.5)
#TOP10
seurat.obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#top2
seurat.obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top2
HP2<-DoHeatmap(seurat.obj.f, features = c(top2$gene)) + NoLegend()
# save. file names
markers.file = paste(save.dir,"/","markers.csv",sep="")
top2.file = paste(save.dir,"/","top2..markers.csv",sep="")
top10.file = paste(save.dir,"/","top10.markers.csv",sep="")
top2.HM.file = paste(save.dir,"/","top2.heatmap.png",sep="")
# write
write.csv(seurat.obj.markers,file=markers.file)
write.csv(top2,file=top2.file)
write.csv(top10,file=top10.file)
ggplot2::ggsave(filename = top2.HM.file, plot = HP2)

###################################################################################
# Get number of cells per cluster and per sample of origin
#############################################################################
cell.count<-table(seurat.obj.f@meta.data[["SCT_snn_res.0.4"]], seurat.obj.f@meta.data[["orig.ident"]])
# Get number of cells per cluster and per sample of origin
clusters<-as.integer(levels(seurat.obj.f@meta.data[["SCT_snn_res.0.4"]]))
cluster.num=max(clusters)+1
cell.count.norm <- matrix(c(rep(0,cluster.num*samples.num)), ncol=samples.num, byrow=TRUE)

colnames(cell.count.norm) <- c("%KOR","%WTR")
rownames(cell.count.norm) <- c(0:(cluster.num-1))
cell.count.norm <- as.table(cell.count.norm)
for (i in 1:samples.num){
  cell.count.norm[,i]<-(cell.count[,i]/sum(cell.count[,i]))*100
}

final_count<-cbind(cell.count,cell.count.norm)
final.count.f<-as.data.frame(final_count)
final.count.f$cluster<-c(0:max(clusters))
cell.count.file = paste(save.dir,"/","cell_count.xlsx",sep="")
rio::export(final.count.f,file=cell.count.file,overwrite=TRUE )
# plot
cell.count.norm<-t(cell.count.norm)
barplot(cell.count.norm,col = c("#1b98e0", "#353436"),beside = TRUE,
        xlab="cluster",ylab="% of cells",main="%cells per cluster",legend.text=TRUE)
#######################################################################
# SAVE data
#######################################################################
saveRDS(new.object,file=seurat_out_file)
sessionInfo() %>% capture.output(file="session_info.txt")
#######################################################################