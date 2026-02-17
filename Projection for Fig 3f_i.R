library(Seurat)
library(scCustomize)
library(ggplot2)
library(dplyr)

# Reading rostarl NMC cells
rostral <- readRDS('rostral.rds')

# Reading human reference and subsetting p.c.w. 5.5 cells
ref <- readRDS('HumanReference.rds')
ref <- SetIdent(ref, value = 'age')
ref <- subset(ref, idents = c('5', '5.5'))
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000)
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
ref <- ScaleData(ref, features = rownames(ref))
ref <- RunPCA(ref, features = VariableFeatures(object = ref))
ref <- RunUMAP(ref, dims = 1:30, reduction = "pca", return.model = TRUE)
ref <- RunUMAP(ref, dims = 1:30, reduction = "pca", return.model = FALSE, reduction.name = "ref.umap")
ref <- SetIdent(ref, value = 'Tissue')

# Fig 3e
DimPlot_scCustom(ref, reduction = 'umap', label = FALSE, colors_use = c("red","#0072B2","#009E73","orange","#009E73"), pt.size =  4) + NoLegend() +NoAxes()+ ggtitle(NULL) # 1750x2000

# Projection
anchors <- FindTransferAnchors(reference = ref, query = rostral, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$Region, dims = 1:30)
rostral <- AddMetaData(rostral, metadata = predictions)
rostral <- MapQuery(anchorset = anchors, reference = ref, query = rostral, reference.reduction = "pca", reduction.model = "umap", new.reduction.name = 'umap')


#######################################################################################
################################ Plotting predicted ID (Fig 3f&g) #####################
#######################################################################################



# Reference Data
ref_umap_coords <- Embeddings(ref, reduction = "umap") %>%
  as.data.frame() %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2) %>% # Standardize column names
  mutate(dataset = "Reference")

# Query Data. Define whether to plot WT or KO
query_seurat_object <- subset(rostral, ident = 'KOR')


query_umap_coords <- Embeddings(query_seurat_object, reduction = "ref.umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  rename(UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2) %>% # Standardize column names
  mutate(dataset = "Query") %>%
  left_join(query_seurat_object@meta.data %>%
              select(predicted.id) %>%
              tibble::rownames_to_column(var = "cell_id"), # Also convert meta.data row names to a column
            by = "cell_id")

# Ensure predicted.id is a factor for consistent coloring
query_umap_coords$predicted.id <- factor(query_umap_coords$predicted.id)

# --- 2. Combine Data Frames ---
combined_umap_data <- bind_rows(ref_umap_coords, query_umap_coords)

# --- 3. Plotting with ggplot2 ---

# Define the order of plotting (Reference first, then Query on top)
combined_umap_data$dataset <- factor(combined_umap_data$dataset, levels = c("Reference", "Query"))

# Get colors for predicted IDs (optional, but good practice for consistency)
# If you have too many predicted IDs, consider using color-blind friendly palettes
num_predicted_ids <- length(levels(query_umap_coords$predicted.id))
colors_for_predicted <- scales::hue_pal()(num_predicted_ids) # Generates distinct colors

# Create the plot
ggplot(combined_umap_data, aes(x = UMAP_1, y = UMAP_2)) +
  # Layer for Reference Data (background, gray)
  geom_point(data = filter(combined_umap_data, dataset == "Reference"),
             color = "gray80", # Light gray for reference background
             size = 3, # Small points for background
             alpha = 0.7) + # Slightly transparent
  
  # Layer for Query Data (on top, colored by predicted identity)
  geom_point(data = filter(combined_umap_data, dataset == "Query"),
             aes(color = predicted.id), # Color by predicted identity
             size = 4, # Slightly larger points for query
             alpha = 0.9) +
  scale_color_manual(values = c("red","orange","#0072B2","#009E73")) + 
  theme_classic() + NoLegend()
# 2000 x 2000



#######################################################################################
################################ Plotting is cluster1 (Fig 3h&i) ######################
#######################################################################################
is.cluster1 <- rostral$seurat_clusters == '1'
rostral <- AddMetaData(rostral, is.cluster1, col.name = 'is.cluster1')

# Reference Data
ref_umap_coords <- Embeddings(ref, reduction = "umap") %>%
  as.data.frame() %>%
  rename(UMAP_1 = umap_1, UMAP_2 = umap_2) %>% # Standardize column names
  mutate(dataset = "Reference")

# Query Data. Define whether to plot WT or KO
query_seurat_object <- subset(rostral, ident = 'WTR')

query_umap_coords <- Embeddings(query_seurat_object, reduction = "ref.umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  rename(UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2) %>% # Standardize column names
  mutate(dataset = "Query") %>%
  left_join(query_seurat_object@meta.data %>%
              select(is.cluster1) %>%
              tibble::rownames_to_column(var = "cell_id"), # Also convert meta.data row names to a column
            by = "cell_id")

# Ensure predicted.id is a factor for consistent coloring
query_umap_coords$is.cluster1 <- factor(query_umap_coords$is.cluster1)

# --- 2. Combine Data Frames ---
combined_umap_data <- bind_rows(ref_umap_coords, query_umap_coords)

# --- 3. Plotting with ggplot2 ---

# Define the order of plotting (Reference first, then Query on top)
combined_umap_data$dataset <- factor(combined_umap_data$dataset, levels = c("Reference", "Query"))
predicted_ids_to_top <- c("FALSE")

combined_umap_data_ordered <- combined_umap_data %>%
  arrange(
    # First, sort by dataset: Reference cells drawn first, then Query cells
    dataset,
    # Then, for Query cells, prioritize specific predicted IDs
    # `predicted.id %in% predicted_ids_to_top` will create TRUE for priority cells, FALSE otherwise.
    # `desc()` makes TRUE (priority cells) come *after* FALSE (other cells)
    # This ensures priority cells are at the end of the query group, thus drawn last.
    desc(is.cluster1 %in% predicted_ids_to_top)
  )



ggplot(combined_umap_data_ordered, aes(x = UMAP_1, y = UMAP_2)) +
  # Layer for Reference Data (background, gray) - filtered here to ensure it's drawn first
  geom_point(data = filter(combined_umap_data_ordered, dataset == "Reference"),
             color = "gray80",
             size = 3,
             alpha = 0.7) +
  
  # Layer for Query Data (on top, colored by predicted identity) - also filtered
  geom_point(data = filter(combined_umap_data_ordered, dataset == "Query"),
             aes(color = is.cluster1),
             size = 4,
             alpha = 0.9) +
  
  # Color scale (choose one of your preferred methods from the previous answer)
  scale_color_manual(values = c('#7DA5B3', "#BE1E2D")) + 
  theme_minimal() + theme_classic() + NoLegend()

# 2000 x 2000


