library(Seurat)
library(scCustomize)
###########################################################
# upload external R file with functions
###########################################################
source("seurat_functions_V2.R")
###########################################################
# convert anndata object to Seurat
###########################################################
sce <- readH5AD("b40faec8-21d7-4d6c-aa69-4146503d3c64.h5ad")
seurat.sub <- as.Seurat(sce,counts = "X", data = "X")
cluster_stats <- Cluster_Stats_All_Samples(seurat_subect = seurat.sub, group_by_var = "tissue")
###########################################################
#continue analysis
###########################################################
Idents(seurat.sub)<-seurat.sub$assay
sub<-subset(seurat.sub,idents="10x 3' v3")
##############################################
# I create a new object with modified gene symbols. Instead of EnsemblIDs
# Get new gene names
new_gene_names <- sub@assays[["originalexp"]]@meta.features[["Gene"]]
counts <- GetAssayData(sub, layer = "counts")
rownames(counts) <- new_gene_names
new_obj <- CreateSeuratObject(counts = counts, meta.data = sub@meta.data)
sub<-new_obj
#############################################################################
#clean
sub<-Add_Cell_QC_Metrics(sub, species = 'human',add_cell_cycle=F)
Idents(sub)<-sub$assay
VlnPlot(sub  , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
        ncol = 3,pt.size=0)
processed_list <- MAD_filtration(sub)
sub.1<-processed_list[[1]]
VlnPlot(sub.1  , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), 
        ncol = 3,pt.size=0)
#
#remove head
Idents(sub.1)<-sub$Region
sub.1<-subset(sub.1,idents = c("Head","Brain"),invert=T)
Idents(sub.1)<-sub.1$Age
levels(sub.1$Age)
sub.1<-subset(sub.1,idents = c("5","12","13","14"),invert=T)
# save 
writeRDS(seurat.sub,"/home/labs/olenderlab/lvzvia/MyRScripts/scRNA_pipline/human_embryonic_atlas/ClassRadialGlia_all_v3_sub.rds")
####################################################################################
## Remove mitochondrial genes and ribosomal genes
sub.1<-remove_ribosomal(sub.1,'human')
#remove ribosomal
mitho_genes <- grep("MT-", rownames(sub.1), value = TRUE)
sub.1 <- subset(sub.1, features = setdiff(rownames(sub.1), mitho_genes))
Idents(sub.1)<-sub$Region
####
# by tissue
sub.1$Subregion <-droplevels(factor(sub.1$Subregion))
table(sub.1$Subregion)
sub.1$tissue <-droplevels(factor(sub.1$tissue))
table(sub.1$tissue)
Idents(sub.1)<-sub.1$tissue

set.seed(123)
clust <- Idents(sub.1)
n_per <- table(clust)

target_total <- 1000 * length(levels(clust))  # same overall budget
min_per_type <- 300    # floor for stability (tune: 100–500)
max_per_type <- 2000   # cap for very large types (tune)

# start with proportional targets
props <- n_per / sum(n_per)
raw <- props * target_total

# apply floor & cap
desired <- pmax(min_per_type, pmin(round(raw), max_per_type))
# adjust to hit exact total
adj_diff <- target_total - sum(desired)

if (adj_diff != 0) {
  # candidates we can still increase/decrease without violating cap/floor
  room_up   <- pmax(0, max_per_type - desired)
  room_down <- pmax(0, desired - min_per_type)
  
  if (adj_diff > 0) {
    # distribute extra to clusters with largest (raw - desired) where room_up > 0
    gain_score <- (raw - desired); gain_score[room_up == 0] <- -Inf
    order_up <- order(gain_score, decreasing = TRUE)
    i <- 1
    while (adj_diff > 0 && i <= length(order_up)) {
      idx <- order_up[i]
      bump <- min(adj_diff, room_up[idx])
      desired[idx] <- desired[idx] + bump
      adj_diff <- adj_diff - bump
      i <- i + 1
    }
  } else {
    # remove from clusters with smallest (raw - desired) where room_down > 0
    loss_score <- (raw - desired); loss_score[room_down == 0] <- Inf
    order_dn <- order(loss_score, decreasing = FALSE)
    i <- 1
    while (adj_diff < 0 && i <= length(order_dn)) {
      idx <- order_dn[i]
      cut <- min(-adj_diff, room_down[idx])
      desired[idx] <- desired[idx] - cut
      adj_diff <- adj_diff + cut
      i <- i + 1
    }
  }
}


desired <- setNames(as.integer(desired), names(n_per))  # reattach names

cells_keep <- unlist(lapply(names(n_per), function(ct) {
  cs <- WhichCells(sub.1, idents = ct)
  k <- min(desired[ct], length(cs))  # use [ ] to keep it simple
  if (k <= 0) character(0) else sample(cs, k)
}), use.names = FALSE)

sub.2 <- subset(sub.1, cells = cells_keep)
sub.2<-NormalizeData(sub.2)
sub.2 <- FindVariableFeatures(sub.2, selection.method = "vst",
                            nfeatures = 20000, verbose = FALSE)
hvgs <- VariableFeatures(sub.2)

counts <- GetAssayData(sub.2 ,  slot = "counts") 
counts<- counts[rownames(counts) %in% hvgs, , drop = FALSE]

# export the filtered raw count matrix (non-normalized)
counts_df <- data.frame(Gene = rownames(counts),
                        as.matrix(counts),
                        check.names = FALSE)
write.table(counts_df,
            file = "reference_miri_scRNA.3_tissue.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# export the labels file (cell_id <tab> cell_type)
sub.2$tissue <- gsub(" ", "_", sub.2$tissue)
levels(factor(sub.2$tissue))
labels_df <- data.frame(cell_id = colnames(sub.2),
                        cell_type = Idents(sub.2))
write.table(labels_df,
            file = "scRNA_labels.3.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
####
# by tissue
labels_df <- data.frame(cell_id = colnames(sub.2),
                        cell_type = sub.2$tissue)
write.table(labels_df,
            file = "scRNA_labels.3_tissue.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

###########################################################
