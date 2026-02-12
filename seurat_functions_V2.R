library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(jsonlite)
library(SingleCellExperiment)
library(BiocSingular)
library(scDblFinder)
library(scater)

#####################################################################################################
# FUNCTIONS
#####################################################################################################
create_folders <-function(save.dir,pca.num,res.num){
  #create main folder
  ifelse(!dir.exists(save.dir), dir.create(save.dir, showWarnings = FALSE), FALSE)
  # cerate folders for the results
  folder.1=paste(save.dir,"/1_QC_plots",sep="")
  ifelse(!dir.exists(folder.1), dir.create(folder.1, showWarnings = FALSE), FALSE)
  #
  solution.folder=paste(save.dir,"/solution_pca",pca.num,"_res",res.num,sep="")
  ifelse(!dir.exists(solution.folder), dir.create(solution.folder, showWarnings = FALSE), FALSE)
  #
  folder.2=paste(solution.folder,"/2_markers",sep="")
  ifelse(!dir.exists(folder.2), dir.create(folder.2, showWarnings = TRUE), FALSE)
  #
  folder.3=paste(solution.folder,"/3_automatic_annotation",sep="")
  ifelse(!dir.exists(folder.3), dir.create(folder.3, showWarnings = FALSE), FALSE)
  #
  folder.4=paste(solution.folder,"/4_comparison_across_conditions",sep="")
  ifelse(!dir.exists(folder.4), dir.create(folder.4, showWarnings = FALSE), FALSE)
  return(solution.folder)
}
dublet_detection_l<-function(the.object){
  the.object<-NormalizeData(the.object)
  sce <- as.SingleCellExperiment(the.object)
  sce <- logNormCounts(sce)
  # Split
  n <- ncol(sce)
  grp <- ceiling(seq_len(n) / 15000)
  sce$chunk <- paste0("chunk_", grp)
  
  
  results <- lapply(unique(sce$chunk), function(ch) {
    tmp <- sce[, sce$chunk == ch]
    set.seed(100 + as.integer(sub("chunk_", "", ch)))
    tmp <- scDblFinder(tmp)
    
    # Clean up potential conflicting metadata
    if ("scDblFinder.selected" %in% names(rowData(tmp))) {
      rowData(tmp)$scDblFinder.selected <- NULL
    }
    
    tmp
  })
  # Combine
  sce <- do.call(cbind, results)
  
  #set.seed(100)
  #sce <- scDblFinder(sce,assays = list(counts = raw_counts))
  
  the.object$doublet_class <- sce$scDblFinder.class
  #the.object <- subset(the.object, subset = doublet_class == "singlet")
  return(the.object)
}
dublet_detection_s<-function(the.object){
  the.object<-NormalizeData(the.object)
  sce <- as.SingleCellExperiment(the.object)
  sce <- logNormCounts(sce)
  sce <- scDblFinder(sce, dbr=0.1)
  the.object$doublet_class <- sce$scDblFinder.class
  #the.object <- subset(the.object, subset = doublet_class == "singlet")
  return(the.object)
}
#####################################################################################################
get_top_PC_genes<-function(the.object){
  pca_results <- the.object[["pca"]]
  loadings <- pca_results@feature.loadings
  top_genes <- loadings[, 1:5] %>% as.data.frame()  # Convert to data frame
  top_genes<-as.data.frame(top_genes)
  top_genes$names=row.names(top_genes)
  # Sort the data frame by the second column in descending order and extract names
  top_n <- 10 # Specify how many top values you want to extract
  top_names<-vector()
  
  for(i in 1:5){
    pc.num=paste("PC_",i,sep="")
    # Sort the data frame by the 'Age' column in descending order
    top_genes <- top_genes %>%
      arrange(desc(!!sym(pc.num)))  # Use desc() for descending order
    top_names<-c(top_names,top_genes$names[1:top_n])
    # Sort the data frame by the 'Age' column in ascending order
    top_genes <- top_genes %>%
      arrange(!!sym(pc.num))  # Use sym and !! to refer to the column
    top_names<-c(top_names,top_genes$names[1:top_n])
  }
  
  return(top_names)
  
}

remove_genes_from_matrix<-function(the.object,genes.file){
  ####################################################################################
  # Remove ribosomal genes
  ###################################################################################
  cc.genes<-readLines(genes.file)
  cc.genes<-as.character(cc.genes)
  all.features = rownames(x = the.object)
  rowNum<-c()
  gene.nums<-as.integer(length(cc.genes))
  endL=gene.nums-1
 
  for(i in 1:gene.nums){
    print(i)
    print(cc.genes[i])
    rowNum[i]<- which(all.features ==cc.genes[i])
    
  }
  all.features <- all.features[-rowNum]
  
  the.object<-subset(x = the.object, features =all.features )
  return(the.object)
}


cluster_cells<-function(the.object,pcas,res,save.dir,to.save,object.reduction,plot.reduction){
  set.seed(1234)
  the.object <- FindNeighbors(the.object, dims = 1:pcas,reduction =object.reduction)
  the.object <- FindClusters(the.object, resolution = res )
  p1<-DimPlot_scCustom(the.object,label.size = 6, color_seed = 2,reduction =plot.reduction)
  p2<-DimPlot_scCustom(the.object,group.by = 'orig.ident',reduction =plot.reduction)
  print (p1|p2)
  if(to.save == 1){
    umap.file=paste(save.dir,"/umap_",pcas,"PCAs_",res,"res.png",sep='')
    ggplot2::ggsave(filename = umap.file, plot = p1|p2,width = 20, height = 15)
  }
  return(the.object)
}
find_markers<-function(the.object,save.dir,normalization.method,pca.num,res.num){
  
  #the.object <- JoinLayers(the.object)
  
  the.object.markers <- FindAllMarkers(the.object, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)#,recorrect_umi=FALSE
  the.object.markers <-subset(the.object.markers,the.object.markers$p_val_adj<0.05)
  #sort
  the.object.markers <- the.object.markers %>%
    arrange(cluster, desc(avg_log2FC))
  #TOP10
  top10 <- the.object.markers %>%
    filter(avg_log2FC > 0) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10)
  #top5
  top5 <- the.object.markers %>%
    filter(avg_log2FC > 0) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5)
  
  HP2<-DoHeatmap(the.object, features = c(top5$gene)) + NoLegend()
  # save. file names
 
  top5.file = paste(save.dir,"/","top5.markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  top10.file = paste(save.dir,"/","top10.markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  top5.HM.file = paste(save.dir,"/","top5.markers.heatmap_",pca.num,"PCA_res_",res.num,".png",sep="")
  
  # write
  markers.file=paste(save.dir,"/","cell_markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  write.csv(the.object.markers,file=markers.file)
  write.csv(top5,file=top5.file)
  write.csv(top10,file=top10.file)
  ggplot2::ggsave(filename = top5.HM.file, plot = HP2,units ="cm",width=30,height=20)
  
}
build.seurat.object<-function(experimentLocation){
  experiment.data <- Read10X(data.dir =experimentLocation)
  seurat.obj<- CreateSeuratObject(counts  = experiment.data, project =projectName , min.cells = 1, min.features = 1)
  
  # read file with mitochondrial genes
  mt.file<-read.csv(file="/home/labs/dvirgur/Collaboration/Medaka/medaka_genome_files/medaka_mitochondrialGenes.csv",header=T)
  mt.genes<-as.vector(mt.file$Gene_name)
  seurat.obj[["percent_mito"]] <- PercentageFeatureSet(seurat.obj , features = mt.genes)
  
  seurat.obj <- Add_Cell_Complexity_Seurat(seurat_object = seurat.obj)
  return(seurat.obj)
}
###################################################################################
# preparing cell count table
####################################################################################
make.countcells.table<-function(the.obj,folder.2,by_phase,test_by){
  if(by_phase==1){
    #table by phase
    cell.count <- Cluster_Stats_All_Samples(seurat_object = the.obj, group_by_var = "Phase")
    cell.count.file = paste(folder.2,"/","cell.count_stage_byPhase_",pca.num,"PCA_res_",res.num,".csv",sep="")
    write.csv(cell.count,file=cell.count.file)
  }
  #table by genotype
  #
  cell.count <- Cluster_Stats_All_Samples(seurat_object = the.obj, group_by_var = test_by)
  cell.count.file = paste(folder.2,"/","cell.count_stage_",test_by,"_",pca.num,"PCA_res_",res.num,".csv",sep="")
  write.csv(cell.count,file=cell.count.file)
  #cell type proportion file
  cell_count_proportions=Proportion_Plot(seurat_object = the.obj, plot_type = "bar",  plot_scale = "percent",split.by ="orig.ident")
  ggsave(paste(folder.2,"/cell_count_proportions.png",sep=""), plot = cell_count_proportions, width = 15, height = 8, dpi = 150)
  # quality parameters per cluster
  qc_per_cluster=VlnPlot_scCustom(the.obj,features = c("nFeature_RNA","percent_oxphos","percent_apop",
                                                          "percent_dna_repair","percent_ieg","percent_hemo"),pt.size = 0.1)
  
  ggsave(paste(folder.2,"/QC_per_cluster.png",sep=""), plot = qc_per_cluster, width = 20, height = 15, dpi = 150)
  
 
}
###################################################################################
# label transfer of annotation using Tabula Moris data set
###################################################################################
run_tabula_moris<-function(query.object,pca.num,normalization_method){
  if(normalization_method == "SCT"){
    ref.obj=readRDS("/home/labs/olenderlab/lvzvia/MyRScripts/scRNA_pipline/TabulaMuris/tabula_muris_facs_12Dec24_SCT.rds")
  }else if(normalization_method == "RNA"){
    ref.obj=readRDS("/home/labs/olenderlab/lvzvia/MyRScripts/scRNA_pipline/TabulaMuris/tabula_muris_facs_12Dec24_RNA.rds")
  }
  Idents(ref.obj)<-ref.obj$cell_ontology_class
 
  if(normalization_method =="SCT"){
    ref.obj<- PrepSCTFindMarkers(ref.obj)
    experiment.anchors <- FindTransferAnchors(reference = ref.obj, query = query.object, dims = 1:pca.num,
                                            reference.reduction = "pca",normalization.method = "SCT")#,recompute.residuals=FALSE)
  }else if(normalization_method =="RNA"){
    experiment.anchors <- FindTransferAnchors(reference = ref.obj, query = query.object, dims = 1:pca.num,
                                              reference.reduction = "pca",normalization.method = "LogNormalize")#,recompute.residuals=FALSE)
  }
  predictions <- TransferData(anchorset = experiment.anchors, refdata = ref.obj@meta.data[["cell_ontology_class"]], dims = 1:15)
  query.object <- AddMetaData(query.object, metadata = predictions)
  Idents(query.object) <- query.object$predicted.id
  DimPlot_scCustom(query.object,pt.size=1,label=F)
  FeaturePlot_scCustom(query.object,features = "prediction.score.max",pt.size=1)
  ########################################################################################################################
  # map query on reference UMAP
  #query.object <- MapQuery(anchorset = experiment.anchors, reference = ref.obj, query = query.object,
  #                         refdata=list(celltype ="cell_ontology_class"),reference.reduction = "pca",
  #                         new.reduction.name="umap.ref",reduction.model = "umap")
  #p1 <- DimPlot(ref.obj, reduction = "umap",  label = TRUE, label.size = 3,
  #              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  #p2 <- DimPlot(query.object, reduction = "umap.ref", group.by = "seurat_clusters", label = TRUE,
  #              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  #print(p1|p2)
  return(query.object)
}
###################################################################################
# Quantile filtration
####################################################################################
quantile_filtration<-function(the.object){
  mt_cutoff <- quantile(the.object@meta.data$percent_mito, 0.90)
  nCount_h_cutoff <- quantile(the.object@meta.data$nCount_RNA, 0.99)
  nCount_l_cutoff <- quantile(the.object@meta.data$nCount_RNA, 0.01)
  nFeature_h_cutoff <- quantile(the.object@meta.data$nFeature_RNA, 0.99)
  nFeature_l_cutoff <- quantile(the.object@meta.data$nFeature_RNA, 0.01)
  cat(mt_cutoff,nFeature_l_cutoff,nFeature_h_cutoff,nCount_l_cutoff,nCount_h_cutoff,"\n")
  final.object <- subset(the.object, subset = percent_mito < mt_cutoff &
                           nCount_RNA < nCount_h_cutoff &
                           nCount_RNA > nCount_l_cutoff &
                           nFeature_RNA < nFeature_h_cutoff &
                           nFeature_RNA > nFeature_l_cutoff)
  return(final.object)
}
###################################################################################
# Custom filtration
###################################################################################
custom_filtration<-function(the.object,mt_percent,Count_genes,Count_UMI){
  
  the.object.f<- subset(the.object, subset = percent_mito < mt_percent &
                          nCount_RNA > Count_UMI[1] & # min UMI
                          nCount_RNA < Count_UMI[2] & # max UMI
                          nFeature_RNA > Count_genes[1] & #min genes
                          nFeature_RNA < Count_genes[2] #max genes
                          )
  return(the.object.f)
}
###################################################################################
# median absolute deviations (MAD) filtration
###################################################################################
MAD_filtration<-function(the.object){
  # this just calculates the parameters.
  # the filtration is performed in the main script
  mad_mt <- mad(the.object$percent_mito)
  median_mt<-median(the.object$percent_mito)
  up_mt=median_mt+(mad_mt*5)
  #
  mad_genes <- mad(the.object$nFeature_RNA)
  median_genes<-median(the.object$nFeature_RNA)
  up_genes=median_genes+(mad_genes*5)
  #
  mad_umi <- mad(the.object$nCount_RNA)
  med_umi=median(the.object$nCount_RNA)
  up_umi=med_umi+(mad_umi*5)
  # filter
  print(up_mt)
  print(up_umi)
  print(up_genes)
  the.object.f<- subset(the.object, subset = percent_mito < up_mt &
                          nCount_RNA > 1000 & # min UMI
                          nCount_RNA < up_umi & # max UMI
                          nFeature_RNA > 500 & #min genes
                          nFeature_RNA < up_genes #max genes
  )

    return(list(the.object.f,up_mt=up_mt,up_umi=up_umi,up_genes=up_genes))
  
}
###################################################################################
# Remove Ribosomal genes
###################################################################################
remove_ribosomal<-function(the.object,organism){
  if(organism == "human"){
    ribo.str=c("^RPS","^RPl")
  }else if(organism == "mouse"){
    ribo.str=c("^Rps","^Rpl")
  }
  
  #list1
  ribo_genes_1 <- grep(ribo.str[1], rownames(the.object), value = TRUE)
  print(ribo_genes_1)
  ngenes=length(ribo_genes_1)
  if(ngenes > 0){
    the.object <- subset(the.object, features = setdiff(rownames(the.object), ribo_genes_1))
  }
  #list2
  ribo_genes_2 <- grep(ribo.str[2], rownames(the.object), value = TRUE)
  print(ribo_genes_2)
  ngenes=length(ribo_genes_1)
  if(ngenes > 0){
    the.object <- subset(the.object, features = setdiff(rownames(the.object), ribo_genes_2))
  }
  return(the.object)
}
###################################################################################
# Do QC plots
###################################################################################
Do_QC_plots<-function(the.object,folder.1,plot_name){
  vln1_1<-VlnPlot(the.object  , features = c("nFeature_RNA"), raster = FALSE,pt.size=0)+ggtitle(paste("nGenes,crude")) + NoLegend()
  vln1_2<-VlnPlot(the.object  , features = c("nCount_RNA"), raster = FALSE,,pt.size=0)+ggtitle(paste("nUMI,crude")) + NoLegend()
  vln1_3<-VlnPlot(the.object  , features = c("percent_mito"), raster = FALSE,pt.size=0)+ggtitle(paste("%mt,crude")) + NoLegend()
  num_cells_report<-as.data.frame(table(the.object$orig.ident))
  text <- tableGrob(num_cells_report)
  grid_plot.1<-grid.arrange(vln1_1,vln1_2,vln1_3,text,ncol = 3)
  #
  file1.name<-paste(folder.1,"/",plot_name,sep="")
  ggsave(file1.name, plot = grid_plot.1, width = 20, height = 30, dpi = 150)
}