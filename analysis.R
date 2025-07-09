################### pre-seting
# install the packages################### pre-seting
# install the packages
install.packages("R.utils")
install.packages("devtools")
devtools::install_github("immunogenomics/presto")
install.packages("BiocManager")
BiocManager::install("glmGamPoi")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scrapper")
BiocManager::install("Nebulosa")
BiocManager::install("pheatmap")
BiocManager::install("msigdbr")
BiocManager::install("clusterProfiler")
BiocManager::install("ComplexHeatmap")
BiocManager::install("circlize")
library(dplyr)
library(Seurat)
library(patchwork)
library(R.utils)
library(Matrix)
library(ggplot2)
library(future)
library(presto)
library(SeuratObject)
library(glmGamPoi)
library(SingleR)
library(celldex)
library(scrapper)
library(Nebulosa)
library(pheatmap)
library(tibble)
library(tidyr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

# set the file location
setwd("file location")

# setting the threads
plan(sequential)

# Set the max-RAM
options(future.globals.maxSize = 200 * 1024^3)  # 200 GiB

################### raw data processing
# Load the PBMC dataset
mat <- readMM("C:/Users/TSH/Desktop/SLE scRNA seq raw data/matrix.mtx")
barcodes <- readLines("C:/Users/TSH/Desktop/SLE scRNA seq raw data/barcodes.tsv")
features <- read.delim("C:/Users/TSH/Desktop/SLE scRNA seq raw data/features.tsv", header = FALSE)
rownames(mat) <- features$V1
colnames(mat) <- barcodes
gene_map <- read.csv("C:/Users/TSH/Desktop/SLE scRNA seq raw data/Gene_Symbol_ID_Convert.csv")
gene_ids <- rownames(mat)
matched_symbols <- gene_map$symbol[match(gene_ids, gene_map$ensembl)]
sum(!is.na(matched_symbols))
rownames(mat) <- ifelse(is.na(matched_symbols), gene_ids, matched_symbols)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = mat, project = "SLE scRNA seq", min.cells = 3, min.features = 200)
metadata <- read.delim("meta.tsv", stringsAsFactors = FALSE)
metadata$disease_state[metadata$disease_state == "na"] <- "healthy"
metadata <- metadata[match(colnames(data), metadata$cellId), ]
data <- AddMetaData(data, metadata = metadata)

# QC processs, calculate the MT genes percentages
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, 
# i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# select high quality identities
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10)

# normalization
data <- NormalizeData(data, verbose = TRUE)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 10000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0 )
plot1 + plot2

# scaling
data <- ScaleData(data)

# PCA
data <- RunPCA(data, features = VariableFeatures(data))

# Examine and visualize PCA results a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca") + NoLegend()
DimHeatmap(data, dims = 1, cells = 5000, balanced = TRUE)
DimHeatmap(data, dims = 1:15, cells = 5000, balanced = TRUE)

# choose the most appropriate dimension of PCA and resolution (0.2~1.5)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 0.5)

# run UMAP
plan("multisession", workers = 12)
data <- RunUMAP(data, dims = 1:15)
DimPlot(data, reduction = "umap")

# find markers for every cluster compared to all remaining cells, report only the positive ones
data.markers <- FindAllMarkers(data, only.pos = TRUE)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(data, features = top10$gene) + NoLegend()

# mapping the each cell clusters
# Blueprint/ENCODE database
BlueprintEncodeData.ref <- celldex::BlueprintEncodeData()
singler.results <- SingleR(test = expr, ref = BlueprintEncodeData.ref, labels = BlueprintEncodeData.ref$label.fine, clusters = clusters)
cell.clusters <- as.character(Idents(data))
cluster.labels <- singler.results$labels
names(cluster.labels) <- rownames(singler.results)
singleR.labels <- cluster.labels[cell.clusters]
names(singleR.labels) <- colnames(data)
data$SingleR.labels.blueprint <- singleR.labels
DimPlot(data, group.by = "SingleR.labels.blueprint", label = TRUE, repel = TRUE) + ggtitle("UMAP annotated by SingleR (Blueprint/ENCODE)") + labs(x = "UMAP1", y = "UMAP2")

################### B cell re-clustering and re-analysis
# select the B cells subset
B_target_cells <- c("naive B-cells", "Class-switched memory B-cells", "Plasma cells")
B_subset_data <- subset(data, subset = SingleR.labels.blueprint %in% B_target_cells)
feature <- c("PHGDH","PSAT1","SHMT1","SHMT2","MTHFD2")
VlnPlot(B_subset_data, features = feature, group.by = "SingleR.labels.blueprint")
DotPlot(B_subset_data, features = feature, group.by = "SingleR.labels.blueprint") + RotatedAxis()
plot_density(B_subset_data, features = feature)

# re-clustering
B_subset_data <- NormalizeData(B_subset_data)
B_subset_data <- FindVariableFeatures(B_subset_data)
B_subset_data <- ScaleData(B_subset_data)
B_subset_data <- RunPCA(B_subset_data)
B_subset_data <- FindNeighbors(B_subset_data, dims = 1:20)
B_subset_data <- FindClusters(B_subset_data, resolution = 0.5)
B_subset_data <- RunUMAP(B_subset_data, dims = 1:15)
DimPlot(B_subset_data, label = TRUE, group.by = "seurat_clusters") + 
  ggtitle("Reclustered B cells")+ labs(x = "UMAP1", y = "UMAP2")

# draw the plot
B_subset_data.markers <- FindAllMarkers(B_subset_data, only.pos = TRUE)
B_subset_data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(B_subset_data, features = top10$gene) + NoLegend()
VlnPlot(B_subset_data, features = feature)
DotPlot(B_subset_data, features = feature) + RotatedAxis()
DotPlot(B_subset_data, features = unique(top10$gene)) + RotatedAxis()
plot_density(B_subset_data, features = feature)


# clusterfying subset data by disease state
flare  <- subset(B_subset_data, subset = disease_state == "flare")
healthy <- subset(B_subset_data, subset = disease_state == "healthy")
managed <- subset(B_subset_data, subset = disease_state == "managed")

# set up a function
plot_gene_density <- function(feature, healthy_obj, flare_obj, managed_obj, mins, maxs, color_palette = c("grey", "gold", "orange", "red", "darkred"),aspect_ratio = 1){
  xlims <- range(
    c(
      Embeddings(healthy_obj, "umap")[,1],
      Embeddings(flare_obj, "umap")[,1],
      Embeddings(managed_obj, "umap")[,1]
    ),
    na.rm = TRUE
  )
  ylims <- range(
    c(
      Embeddings(healthy_obj, "umap")[,2],
      Embeddings(flare_obj, "umap")[,2],
      Embeddings(managed_obj, "umap")[,2]
    ),
    na.rm = TRUE
  )
  color_scale <- scale_color_gradientn(colors = color_palette, limits = c(mins, maxs), oob = scales::squish)
  p1 <- plot_density(healthy_obj, features = feature) + 
    ggtitle(paste0("Healthy - ", feature)) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    coord_cartesian(xlim = xlims, ylim = ylims) +
    theme(aspect.ratio = aspect_ratio)+
    color_scale
  p2 <- plot_density(flare_obj, features = feature) + 
    ggtitle(paste0("Flare - ", feature)) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    coord_cartesian(xlim = xlims, ylim = ylims) +
    theme(aspect.ratio = aspect_ratio)+
    color_scale
  p3 <- plot_density(managed_obj, features = feature) + 
    ggtitle(paste0("Managed - ", feature)) + 
    labs(x = "UMAP1", y = "UMAP2") + 
    coord_cartesian(xlim = xlims, ylim = ylims) +
    theme(aspect.ratio = aspect_ratio)+
    color_scale
  
  return(p1 + p2 + p3)}

# print and save the density plot of genes
p <- plot_gene_density("PHGDH", healthy, flare, managed, 0, 0.3)
print(p)

# re-mapping the clusters
# Blueprint/ENCODE
expr <- GetAssayData(B_subset_data, slot = "data")
clusters <- Idents(B_subset_data)
singler.results <- SingleR(test = expr, ref = BlueprintEncodeData.ref, labels = BlueprintEncodeData.ref$label.fine, clusters = clusters)
cell.clusters <- as.character(Idents(B_subset_data))
cluster.labels <- singler.results$labels
names(cluster.labels) <- rownames(singler.results)
singleR.labels <- cluster.labels[cell.clusters]
names(singleR.labels) <- colnames(B_subset_data)
B_subset_data$SingleR.labels.blueprint <- singleR.labels
DimPlot(B_subset_data, group.by = "SingleR.labels.blueprint", label = TRUE, repel = TRUE) + ggtitle("UMAP annotated by SingleR (Blueprint/ENCODE)") + labs(x = "UMAP1", y = "UMAP2")

# draw the specific gene expression
feature <- c("CD27","CD19","CD38","PHGDH","PSAT1","SHMT1","SHMT2","MTHFD2")
DotPlot(B_subset_data, features = feature, group.by = "seurat_clusters") + RotatedAxis()
DotPlot(B_subset_data, features = unique(top10$gene)) + RotatedAxis()
top10 <- B_subset_data.markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(!grepl("^ENSG", gene)) %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()
DotPlot(B_subset_data, features = unique(top10$gene)) + RotatedAxis()

################### GSEA
# prepare the gene expression matrix
Idents(B_subset_data) <- "disease_state"
deg <- FindMarkers(
  B_subset_data,
  ident.1 = "flare",
  ident.2 = "healthy",
  logfc.threshold = 0,
  min.pct = 0.1,
  only.pos = FALSE)
deg$gene <- rownames(deg)
gene_list <- deg %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(avg_log2FC)
names(gene_list) <- deg %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

# prepare gsea database (Hallmark)
geneset_df <- msigdbr(species = "Homo sapiens", category = "H")
term2gene <- geneset_df %>% dplyr::select(gs_name, gene_symbol)

# setting the seed
set.seed(123)

# run GSEA
gsea_result <- GSEA(geneList = gene_list, TERM2GENE = term2gene, pvalueCutoff = 0.05)

# draw the plot
dotplot(gsea_result, showCategory = 15)

# save the gsea data
write.csv(gsea_result@result, "gsea_flare_vs_healthy_hallmark.csv")

# prepare gsea database (c5)
msig_c5_gobp <- msigdbr( species = "Homo sapiens", category = "C5", subcollection = "GO:BP")
term2gene_c5 <- msig_c5_gobp %>% dplyr::select(gs_name, gene_symbol)

# run GSEA
gsea_result <- GSEA(geneList = gene_list, TERM2GENE = term2gene_c5, pvalueCutoff = 0.05, eps = 0)

# draw the plot
dotplot(gsea_result, showCategory = 15)

# save the gsea data
write.csv(gsea_result@result, "gsea_flare_vs_healthy_c5.csv")

# volcano plot
deg$significance <- "Not Sig"
deg$significance[deg$p_val_adj < 0.05 & deg$avg_log2FC > 1] <- "Up"
deg$significance[deg$p_val_adj < 0.05 & deg$avg_log2FC < -1] <- "Down"
deg <- deg %>%
  mutate(log10_pval = -log10(pmax(p_val_adj, 1e-300)))
top_genes <- deg %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(deg, aes(x = avg_log2FC, y = log10_pval, color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 100) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  ylim(0, 350)+
  theme_minimal()+
  labs(x = "log2 Fold Change",
       y = "-log10 adjusted p-value",
       title = "Volcano Plot: flare vs healthy")

###################
# select the plasma cells subset
plasma_cell_subset_data <- subset(B_subset_data, subset = SingleR.labels.blueprint == "Plasma cells")

# clusterfying subset data by disease state
flare  <- subset(plasma_cell_subset_data, subset = disease_state == "flare")
healthy <- subset(plasma_cell_subset_data, subset = disease_state == "healthy")
managed <- subset(plasma_cell_subset_data, subset = disease_state == "managed")

# prepare the gene expression matrix
Idents(plasma_cell_subset_data) <- "disease_state"
deg <- FindMarkers(
  plasma_cell_subset_data,
  ident.1 = "flare",
  ident.2 = "healthy",
  logfc.threshold = 0,
  min.pct = 0.1,
  only.pos = FALSE)
deg$gene <- rownames(deg)
gene_list <- deg %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(avg_log2FC)
names(gene_list) <- deg %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  pull(gene)

# run GSEA
gsea_result <- GSEA(geneList = gene_list, TERM2GENE = term2gene, pvalueCutoff = 0.05)

# draw the plot
dotplot(gsea_result, showCategory = 15)

# run GSEA
gsea_result <- GSEA(geneList = gene_list, TERM2GENE = term2gene_c5, pvalueCutoff = 0.05, eps = 0)

# draw the plot
dotplot(gsea_result, showCategory = 15)

# volcano plot
deg$significance <- "Not Sig"
deg$significance[deg$p_val_adj < 0.05 & deg$avg_log2FC > 1] <- "Up"
deg$significance[deg$p_val_adj < 0.05 & deg$avg_log2FC < -1] <- "Down"
top_genes <- deg %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(20)
ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 100) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  ylim(0, 30)+
  theme_minimal()+
  labs(x = "log2 Fold Change",
       y = "-log10 adjusted p-value",
       title = "Volcano Plot: flare vs healthy")