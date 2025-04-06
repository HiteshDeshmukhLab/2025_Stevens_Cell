# Code used to analyze and present the human scRNAseq data included in Stevens et al., Cell, 2025 

rm(list=ls())
library(tidyverse)
library(seuratTools)
library(Seurat)
library(Matrix)
library(DropletUtils)
library(reshape2)
library(SeuratWrappers)

# Define the directory containing subfolders for each sample
data_directory <- "/path"
setwd("/path")

# Get a list of subdirectories
sample_dirs <- list.dirs(data_directory, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store Seurat objects
seurat_list <- list()

# Loop through each subdirectory and create Seurat objects

for (sample_dir in sample_dirs) {
  # Extract the folder name to use as the Seurat object name
  sample_name <- basename(sample_dir)
  # Read the 10x data from the current folder
  data <- Read10X(data.dir = sample_dir, gene.column = 1)
  # Create a Seurat object using the sample name as the project name
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  # Add sample name as metadata
  seurat_obj$sample <- sample_name
  # Add `exp` column based on folder name
  if (grepl("Control", sample_name)) {
    seurat_obj$exp <- "Control"
  } else if (grepl("Dysbiosis", sample_name)) {
    seurat_obj$exp <- "Dysbiosis"
  } else {
    seurat_obj$exp <- "Unknown"  # Optional: Label for unmatched cases
  }
  # Assign the Seurat object to the list using the folder name
  seurat_list[[sample_name]] <- seurat_obj
}

# Merge all Seurat objects into one combined object
combined_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list)
)

# Annotate MT, RPS, and RPL genes
seu [['percent.mt']] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu [['percent.rps']] <- PercentageFeatureSet(seu, pattern = "^RPS")
seu [['percent.rpl']] <- PercentageFeatureSet(seu, pattern = "^RPL")
seu$percent.rp <- seu$percent.rps + seu$percent.rpl
seu # 37182 features across 30516 samples within 1 assay
head(seu[[]])

# Filter low quality cells.
minFeature <- 200
maxFeature <- 7500
minCount <- 400
maxCount <- 40000
maxMT <- 5
seu <- subset(seu, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount &
                nCount_RNA < maxCount & percent.mt < maxMT)
VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt"))

# Filter out additional genes: filter genes requiring a min.value (log-normalized) in at least min.cells
# Expression of 1 in at least 200 cells.
FilterGenes <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    genes.use <- rownames(object)
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- GetAssayData(object)[genes.use, ]
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(GetAssayData(object) > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      object = object[genes.use, ]
    }
    object <- LogSeuratCommand(object = object)
    return(object)
  }
seu <- NormalizeData(seu)
seu <- FilterGenes(object = seu, min.value = 1, min.cells = 200) # 400
seu # 10880 features across 25599 samples within 1 assay  

# Calculate Cell-Cycle with Seurat, the list of genes comes with Seurat (only for human)
s.genes <- (cc.genes$s.genes)
g2m.genes <- (cc.genes$g2m.genes)
Genes <-rownames(seu)
s.genes  <- s.genes[s.genes  %in% Genes]
g2m.genes <- g2m.genes[g2m.genes %in% Genes]
seu <- CellCycleScoring(seu, s.features = s.genes,g2m.features = g2m.genes,set.ident = TRUE)

# Prepare object.list for integration
unique(seu$replicate)
object.list <- SplitObject(seu, split.by = 'replicate')
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = object.list)
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors using references (RPCA)
anchors <- FindIntegrationAnchors(object.list = object.list, reduction = 'rpca', dims = 1:20, k.anchor = 30)
seu <- IntegrateData(anchorset = anchors, dims = 1:20)
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:5)
seu <- FindClusters(seu, resolution = 0.1)
seu <- RunUMAP(object = seu, dims = 1:5)
Idents(seu) <- 'seurat_clusters'
DimPlot(seu, label = T, repel = T, split.by = 'exp')+NoLegend()
DimPlot(seu, label = T, repel = T)+NoLegend()
table(seu$replicate, seu$seurat_clusters)
save(seu, file = '/path/seu.all.HUMAN.R1.RDS')
# Eliminate cluster 3, 5 and 6 as they have low number of cells
seu <- subset(seu, idents = c('3','5', '6'), invert = T)

# Rerun standard workflow 
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.15)
seu <- RunUMAP(object = seu, dims = 1:10)
Idents(seu) <- 'seurat_clusters'
DimPlot(seu, label = T, repel = T, split.by = 'exp')+NoLegend()
DimPlot(seu, label = T, repel = T)+NoLegend()
table(seu$replicate, seu$seurat_clusters)

# Save object
seu <- readRDS(file = 'seu.RDS')
DefaultAssay(seu) <- "RNA"
FeaturePlot(seu, features = c('SELL', 'KLF2', 'LEF1', 'CD4', 'CD8A', 'CXCR6','ITGA1', 'GZMB'), label = T)
FeaturePlot(seu, features = c('GZMB', 'IFNG', 'CCL4', 'CCL3', 'NKG7', 'CD8A', 'GNLY', 'KLRD1', 'IL7R', 'VIM', 'CAPG'), label = T) #0 and 2 looks like CD8
FeaturePlot(seu, features = c('TCF7', 'IL4R', 'SELL', 'CD4', 'ID3', 'CCR7'), label = T) # 2 looks like CD4
FeaturePlot(seu, features = c('FOXP3', 'BGS1', 'STAT1', 'LGALS3', 'IL2RA', 'CTLA4','TIGIT', 'IL7R'), label = T) # 5 looks like FOXP3 cells
FeaturePlot(seu, features = c('KLF6', 'ITGA1', 'VIM', 'ITGA1', 'ITGAE', 'CEBPD', 'CAPG'), label = T)
FeaturePlot(seu, features = c('CD4', 'CD8A', 'GZMB', 'GZMK', 'TBX21', 'GATA3', 'PRF1'), label = T)
FeaturePlot(seu, features = c('PDCD1', 'TIGIT', 'CTLA4','LAG3', 'HAVCR2'), label = T) #Exhaustion Signature 0 and 4
FeaturePlot(seu, features = c('CD3D','CD8A', 'HLA-DRA', 'IL17'), label = T)
FeaturePlot(seu, features = c('CD55','HLA-DRB1', 'ACTB', 'HLA-DPB1','SOX4', 'IL6ST', 'NFKBIZ', 'HLA-DRA', 'RUNX1', 'HLA-DPA1', 'HES'),label = T) # Activated T cells
FeaturePlot(seu, features = c('NFIL3', 'TCF7', 'ARNTL', 'CLOCK', 'LEF1'), label = T, min.cutoff = "q01", max.cutoff = "q20") # 1 and 3 are naive T cells

# From Hubmap and # Tissue residency signature from Kumar et al PMID: 28930685 to create Module scores
CD8.MEM <- c('CD8B','CD8A','GZMK','IL7R','BCL11B','CD3E','PIK3IP1','CD6','TRAC','SPOCK2')
CD8.Naive <- c('NELL2','CD8B','S100B','TRABD2A','CCR7','OXNAD1','FLT3LG','MAL','LEF1','CAMK4')
CD4 <- c('MAL','LEF1','IL7R','TCF7','LTB','LEPROTL1','CD3E','CD3D','CD3G','CD27')
CD8 <- c('CD8B','CD8A','NELL2','CD3D','CD3E','S100B','GZMH','CD3G','TRGC2','CCL5')
CD4.Effector <- c('TNFRSF4','AQP3','TNFRSF25','MAL','IL7R','TRAT1','RORA','FLT3LG','KLRB1','CD6')
CD8.Effector <- c('GZMK','CD160','CD8A','KLRB1','LAG3','CCL4','XCL2','GZMM','XCL1','CCL3', 'IFNG','GZMB','TIGIT')
CD4.MEM <- c('MAL','TRAT1','AP3M2','CAMK4','AQP3','IL7R','RCAN3','BCL11B','LEF1','FLT3LG')
NK <- c('GNLY', 'TYROBP', 'NKG7', 'FCER1G', 'GZMB', 'TRDC', 'PRF1', 'FGFBP2', 'SPON2', 'KLRF1')
CD4.ACT <- c('IL2', 'ODC1', 'PSAT1', 'WARS', 'PYCR1', 'TNF', 'TCF7', 'IL4R', 'SELL')
TREG <- c('RTKN2', 'FOXP3', 'AC133644.2', 'CD4', 'IL2RA', 'TIGIT', 'CTLA4', 'FCRL3', 'LAIR2', 'IKZF2')
Tissue.Residency <- c ('FAM129A', 'ITGA1', 'KLF6', 'NR3C1M', 'LGALS3', 'CD82', 'AMICA1',
                       'VIM', 'CAPG', 'ITM2A', 'ITM2C', 'TMEM14C', 'PL2', 'PHLDA1', 'ITGAE', 'CEBPD', 'ANKRD28')
FeaturePlot(seu, features = CD8.Naive,label = T)
FeaturePlot(seu, features = CD4.Effector,label = T)
FeaturePlot(seu, features = Tissue.Residency, label = T)
FeaturePlot(seu, features = CD8.Effector, label = T)
FeaturePlot(seu, features = CD4.MEM, label = T)
FeaturePlot(seu, features = c('CCR7', 'TCF7', 'SELL'),label = T) # 0 = CD4 naive
FeaturePlot(seu, features = c('ITGA1', 'GZMM', 'PRF1', 'NKG7'),label = T) # 1 = CD8 TRM
FeaturePlot(seu, features = c('LAG3', 'CXCR6', 'PRF1', 'IL2RA', 'PYCR1'),label = T) # 2 = CD4 effector 
FeaturePlot(seu, features = c('ITGA1', 'GZMM', 'PRF1', 'NKG7', 'LAG3', 'CXCR6', 'GZMB'),label = T) # 4 = CD8 effector

Genes <-rownames(seu)
View(as.matrix(Genes))
Tissue.Residency  <- Tissue.Residency[Tissue.Residency  %in% Genes]
CD8.MEM <- CD8.MEM[CD8.MEM %in% Genes]
CD4.MEM <- CD4.MEM[CD4.MEM %in% Genes]
CD4.Effector <- CD4.Effector[CD4.Effector %in% Genes]
CD8.Effector <- CD8.Effector[CD8.Effector %in% Genes]
CD4.ACT <- CD4.ACT[CD4.ACT %in% Genes]
CD4 <- CD4[CD4 %in% Genes]
CD8 <- CD8[CD8 %in% Genes]
TREG <- TREG[TREG %in% Genes]

DefaultAssay(seu) <-"RNA"
seu <- JoinLayers(seu)
seu <- AddModuleScore(seu, features = list(Tissue.Residency), ctrl.size = 5, name ='TR')
seu <- AddModuleScore(seu, features = list(CD8.MEM), ctrl.size = 5, name ='CD8.MEM')
seu <- AddModuleScore(seu, features = list(CD4.MEM), ctrl.size = 5, name ='CD4.MEM')
seu <- AddModuleScore(seu, features = list(CD4.Effector), ctrl.size = 5, name ='CD4.EFF')
seu <- AddModuleScore(seu, features = list(CD8.Effector), ctrl.size = 5, name ='CD8.EFF')
seu <- AddModuleScore(seu, features = list(CD4.ACT), ctrl.size = 5, name ='CD4.ACT')
seu <- AddModuleScore(seu, features = list(CD4), ctrl.size = 5, name ='CD4')
seu <- AddModuleScore(seu, features = list(CD8), ctrl.size = 5, name ='CD8')
seu <- AddModuleScore(seu, features = list(TREG), ctrl.size = 5, name ='TREG')

# Visualized metadata slots
head(seu[[]])
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD8.MEM1')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD4.MEM1')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD4.EFF1')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD8.EFF1')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD4.ACT1')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD41')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'CD81')))
seu[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seu, vars = 'TREG1')))

# Plot 
FeaturePlot(object = seu, features = 'TR1', min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, cols = c("grey", "#440154FF"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD8.MEM1', min.cutoff = "q10", max.cutoff = "q90", pt.size = 1, cols = c("grey", "#8c564b"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD4.MEM1', min.cutoff = "q20", max.cutoff = "q90", pt.size = 1, cols = c("grey", "#aa40fc"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD4.EFF1', min.cutoff = "q10", max.cutoff = "q80", pt.size = 1, cols = c("grey", "#ff7f0e"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD8.EFF1', min.cutoff = "q40", max.cutoff = "q80", pt.size = 1, cols = c("grey", "#279e68"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD4.ACT1', min.cutoff = "q20", max.cutoff = "q80", pt.size = 1, cols = c("grey", "#440154FF"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD41', min.cutoff = "q05", max.cutoff = "q60", pt.size = 1, cols = c("grey", "#aa40fc"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'CD81', min.cutoff = "q10", max.cutoff = "q80", pt.size = 1, cols = c("grey", "#aa40fc"), label = T) + NoLegend()
FeaturePlot(object = seu, features = 'TREG1', min.cutoff = "q30", max.cutoff = "q90", pt.size = 1, cols = c("grey", "#440154FF"), label = T) + NoLegend()

# 4 - CD8 Effector, 1 - TRM, 0 - CD4 Mem/Naive, 2 - CD4 Effector\
# Assign cluster IDs 
seu$main.cell <- seu$seurat_clusters
seu$main.cell <- ifelse(seu$seurat_clusters %in% c('1'), 'CD8 TRM', seu$main.cell)
seu$main.cell <- ifelse(seu$seurat_clusters %in% c('0'), 'CD4 Naive', seu$main.cell)
seu$main.cell <- ifelse(seu$seurat_clusters %in% c('2'), 'CD4 Effector', seu$main.cell)
seu$main.cell <- ifelse(seu$seurat_clusters %in% c('4'), 'CD8 Effector', seu$main.cell)
unique(seu$main.cell)
Idents(seu) <- 'main.cell'
DimPlot(seu, label = T)

# Plot cells in each cluster 
# Create a contingency table for clusters and replicates
cluster_rep_table <- table(seu$main.cell, seu$exp, seu$replicate)
# Convert the table to a data frame for easier manipulation
cluster_rep_df <- as.data.frame(cluster_rep_table)
# Rename columns for clarity
colnames(cluster_rep_df) <- c("Cluster", "Exp", "Replicate", "CellCount")
head(cluster_rep_df)
# Calculate the total number of cells for each replicate
total_cells_per_replicate <- cluster_rep_df %>%
  group_by(Replicate) %>%
  summarize(TotalCells = sum(CellCount))
# Merge total cells back into the original data frame
cluster_rep_df <- cluster_rep_df %>%
  left_join(total_cells_per_replicate, by = "Replicate")
# Calculate the percentage of cells in each cluster for each replicate
cluster_rep_df <- cluster_rep_df %>%
  mutate(Percentage = (CellCount / TotalCells) * 100)
# Calculate mean and SEM for each cluster and condition
cluster_stats <- cluster_rep_df %>%
  group_by(Cluster, Exp) %>%
  summarize(
    MeanPercentage = mean(Percentage),
    SEM = sd(Percentage) / sqrt(n()),
    .groups = "drop"
  )

# Create scatter plot with mean and SEM
ggplot() +
  # Plot individual data points
  geom_jitter(data = cluster_rep_df, 
              aes(x = Cluster, y = Percentage, color = Exp), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              alpha = 0.7, size = 2) +
  # Plot mean percentages
  geom_point(data = cluster_stats, 
             aes(x = Cluster, y = MeanPercentage, color = Exp), 
             position = position_dodge(width = 0.5), 
             size = 4, shape = 18) +
  # Plot error bars for SEM
  geom_errorbar(data = cluster_stats, 
                aes(x = Cluster, ymin = MeanPercentage - SEM, ymax = MeanPercentage + SEM, color = Exp), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  labs(title = "Percentage of Cells in Each Cluster by Condition",
       x = "Cluster",
       y = "Percentage (%)",
       color = "Exp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Use Monocle to generate pretty UMAP 
library(SeuratWrappers)
library(monocle3)
DefaultAssay(seu) <- "RNA"
cds <- as.cell_data_set(seu)
cds <- cluster_cells(cds)
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "main.cell", cell_size = 0.5, label_cell_groups = FALSE, alpha = 0.25, cell_stroke = 0.5)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = T)

# find markers for every cluster compared to all remaining cells, report only the positive ones #
Idents(seu) <- 'main.cell'
DefaultAssay(seu) <- 'RNA'
Immune.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(Immune.markers)
Cluster.top.markers <-Immune.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% print(n=250)
View(Cluster.top.markers)
DefaultAssay(seu) <- "RNA"
seu.scale <-ScaleData(seu, features = Cluster.top.markers$gene)
# Remove duplicates for making heatmaps
List <- Cluster.top.markers[!duplicated(Cluster.top.markers$gene), ]
List <- column_to_rownames(List, var = "gene")
List <- rownames(List)
library(viridis)
DoHeatmap(seu.scale, features = List, size = 4) + scale_fill_viridis_b() + NoLegend()

# Create topGOdata object
library(topGO)
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)

# Identify DEG genes between 1) Control and 2) Dysbiosis groups.
Idents(seu) <- "exp"
seu.exp.marker <- FindMarkers(seu, ident.1 = "Control", ident.2 = "Dysbiosis", min.pct = 0.10, logfc.threshold = 1, assay ="RNA")
# view results
View(seu.exp.marker)
genes <- rownames(seu.exp.marker)
gene_names <- genes[, 1] 
gene_vector <- rep(1, length(gene_names))
names(gene_vector) <- gene_names

# Filter results for top 200 DEG
seu.exp.top.markers <- rbind((seu.exp.marker  %>% top_n(n = 100, wt = avg_log2FC)), (seu.exp.marker  %>% top_n(n = -100, wt = avg_log2FC)))
seu.exp.top.markers.NABX <- rbind((seu.exp.marker  %>% top_n(n = 100, wt = avg_log2FC)), (seu.exp.marker  %>% top_n(n = 100, wt = avg_log2FC)))
seu.exp.top.markers.ABX <- rbind((seu.exp.marker  %>% top_n(n = 100, wt = avg_log2FC)), (seu.exp.marker  %>% top_n(n = -100, wt = avg_log2FC)))
View(seu.exp.top.markers)
write.table(seu.exp.top.markers, file = "/path/seu.exp.tsv", sep="\t")

# Save File
library(SeuratDisk)
SaveH5Seurat(seu, filename = "seurat_object.h5Seurat")
Convert("seurat_object.h5Seurat", dest = "h5ad")

# GO Terms enriched in Control vs Dysbiosis

# Create topGOdata object
library(topGO)
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = gene_vector,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)


# Step 1: Define cluster metadata column
cluster_col <- "main.cell"  # Adjust based on your Seurat metadata

# Step 2: Aggregate Counts by Cluster
# Extract raw counts from the Seurat object
counts <- GetAssayData(seurat_object, assay = "RNA", slot = "counts") %>%
  as.matrix()

# Add cluster information to cells
cell_metadata <- seurat_object@meta.data[, cluster_col, drop = FALSE]
cell_metadata$cell <- rownames(cell_metadata)

# Transpose counts to match cell metadata (cells as rows)
counts <- t(counts)
counts <- as.data.frame(counts)
counts$cell <- rownames(counts)

# Merge counts with metadata
counts_with_metadata <- merge(counts, cell_metadata, by = "cell")

# Aggregate counts by cluster
pseudobulk_counts <- counts_with_metadata %>%
  pivot_longer(-c(cell, main.cell), names_to = "gene", values_to = "count") %>%
  group_by(main.cell, gene) %>%
  summarise(count = sum(count), .groups = "drop")

# Step 3: Prepare Data for edgeR
# Pivot pseudobulk counts into a wide matrix
count_matrix <- pseudobulk_counts %>%
  pivot_wider(names_from = main.cell, values_from = count, values_fill = 0) %>%
  column_to_rownames("gene")

# Step 4: Run edgeR for Differential Expression
# Perform exact tests for each cluster
results_list <- lapply(colnames(count_matrix), function(cluster) {
  # Create a group variable: cluster vs others
  group <- ifelse(colnames(count_matrix) == cluster, cluster, "other")
  
  # Create DGEList object
  y <- DGEList(counts = count_matrix, group = group)
  y <- calcNormFactors(y)
  
  # Estimate dispersion
  y <- estimateCommonDisp(y)
  
  # Perform exact test
  res <- exactTest(y, pair = c("other", cluster))
  
  # Extract results
  res_table <- topTags(res, n = Inf)$table
  res_table$gene <- rownames(res_table)
  res_table$cluster <- cluster
  return(res_table)
})

# Combine results from all clusters
all_markers <- do.call(rbind, results_list)

# Step 5: Filter Significant Marker Genes
# Filter for significant markers
significant_markers <- all_markers %>%
  filter(FDR < 0.05) %>%
  arrange(cluster, PValue)

# Save significant markers to a CSV file
write.csv(significant_markers, "/path/significant_markers.csv", row.names = FALSE)
View(significant_markers)

# Step 6: Visualize with a Heatmap
# Extract normalized counts
y <- DGEList(counts = count_matrix)  # Recreate DGEList for normalization
y <- calcNormFactors(y)
normalized_counts <- cpm(y, log = TRUE)

# Filter for significant genes
heatmap_data <- normalized_counts[rownames(normalized_counts) %in% significant_markers$gene, ]

# Scale data
scaled_data <- t(scale(t(heatmap_data)))

# Plot heatmap
pheatmap(
  scaled_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)

# Step 7: Save Heatmap Data (Optional)
write.csv(scaled_data, "/path/heatmap_data.csv", row.names = TRUE)

# Define the marker genes for each subset
highlight_genes <- c(
  # CD4 Naive
  "CCR7", "LEF1", "TCF7", "IL7R", "FOXP1", "IKZF1", "CD27", 
  # CD4 Effector
  "IFNG", "IL13", "IL17F", "TNF", "GZMB", "PRF1", "TBX21", "RORC", "CXCR3",
  # CD8 Naive
  "CCR7", "LEF1", "TCF7", "CD27", "IL7R", "CD28", "FOXP1", "LCK", "CD8A",
  # CD8 Effector
  "GZMB", "PRF1", "IFNG", "TNF", "CXCR3", "KLRG1", "EOMES", "TBX21", "FASLG", "GNLY", "CD8A",
  # CD8 TRM,
  "CD69", "ITGAE", "CXCR6", "ITGB1", "PDCD1", "ZBTB46", "PRDM1", "RUNX3", "GZMK", "CXCL13"
)

# Filter normalized counts for the highlighted genes
highlight_heatmap_data <- normalized_counts[rownames(normalized_counts) %in% highlight_genes, ]

# Check if any genes are missing and warn
missing_genes <- setdiff(highlight_genes, rownames(highlight_heatmap_data))
if (length(missing_genes) > 0) {
  message("The following genes were not found in the dataset: ", paste(missing_genes, collapse = ", "))
}

# Scale the data for visualization
highlight_scaled_data <- t(scale(t(highlight_heatmap_data)))

# Plot the heatmap
pheatmap(
  highlight_scaled_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap of Selected Marker Genes"
)


# Save heatmap data (Optional)
write.csv(highlight_scaled_data, "highlight_heatmap_data.csv", row.names = TRUE)

# Idenify change in propotion of cells # 
# Step 1: Install and Load Required Packages
if (!requireNamespace("scProportionTest", quietly = TRUE)) {
  install.packages("scProportionTest")
}
library(scProportionTest)
library(Seurat)

seu.T <- seurat_object

# Step 2: Verify Metadata in Seurat Object
# Ensure your Seurat object `seu.T` has cluster information and experimental group columns
head(seu.T@meta.data)

# If necessary, rename metadata columns
seu.T@meta.data$custom_clusters <- seu.T@meta.data$main.cell # Cluster labels
seu.T@meta.data$orig.ident <- seu.T@meta.data$replicate # Sample identifiers
seu.T@meta.data$exp <- seu.T@meta.data$exp # Experimental groups

# Step 2: Run Permutation Test
prop_test <- sc_utils(seu.T)

# Step 3: Verify Metadata
# Check the metadata to ensure it contains the required columns
head(prop_test)

# Ensure the following columns exist:
# - "custom_clusters": Cluster identity.
# - "exp": Experimental group (e.g., NABX and ABX).
# - "orig.ident": Sample identifiers.

# If needed, rename columns in the metadata
seu.T@meta.data@meta_data$custom_clusters <- seu.T@meta.data@meta_data$main.cell
seu.T@meta.data@meta_data$exp <- seu.T@meta.data@meta_data$exp
seu.T@meta.data@meta_data$orig.ident <- seu.T@meta.data@meta_data$replicate

# Step 4: Perform Permutation Test
# Run permutation test to compare proportions between NABX and ABX
prop_test <- permutation_test(
  sc_utils_obj = prop_test, # The sc_utils object
  cluster_identity = "custom_clusters", # Column containing cluster labels
  sample_1 = "Dysbiosis", # Group 1 to compare
  sample_2 = "Control", # Group 2 to compare
  sample_identity = "exp" # Column identifying sample replicates
)

# Step 5: View Results
# Print results of the test
permutation_plot(prop_test)

# Extract permutation results
permutation_results <- prop_test@results$permutation
write.csv(permutation_results, "/path/permutation_result.csv", row.names = TRUE)

# Filter significant clusters (pval < 0.05)
significant_clusters <- permutation_results[permutation_results$pval < 0.05, ]

# Print significant clusters
cat("Significant Clusters (pval < 0.05):\n")
print(significant_clusters)

# Extract proportions for Dysbiosis and Control
proportions <- data.frame(
  Cluster = permutation_results$clusters,
  Dysbiosis = permutation_results$Dysbiosis,
  Control = permutation_results$Control
)

# Ensure correct alignment of names and matrix
height_matrix <- as.matrix(proportions[, c("NABX", "ABX")])
cluster_names <- proportions$Cluster  # Cluster names

# Check alignment
if (nrow(height_matrix) != length(cluster_names)) {
  stop("Mismatch between number of rows in height matrix and cluster names.")
}

# Generate the bar plot
# Ensure correct alignment of data
height_matrix <- as.matrix(proportions[, c("NABX", "ABX")])
cluster_names <- proportions$Cluster
p_values <- permutation_results$pval  # Extract p-values

# Check alignment
if (nrow(height_matrix) != length(cluster_names) || nrow(height_matrix) != length(p_values)) {
  stop("Mismatch between data dimensions.")
}

# Generate the bar plot
bar_positions <- barplot(
  height = t(height_matrix),  # Transpose for plotting
  beside = TRUE,
  names.arg = cluster_names,
  col = c("blue", "red"),
  legend.text = c("NABX", "ABX"),
  args.legend = list(x = "topright"),
  main = "Cluster Proportions: NABX vs ABX (with p-values)",
  xlab = "Clusters",
  ylab = "Proportion"
)

# Add p-values above the bars
text(
  x = colMeans(bar_positions),  # Position above bars
  y = apply(height_matrix, 1, max) + 0.02,  # Slightly above the tallest bar
  labels = paste0("p=", signif(permutation_results$pval, digits = 3)),  # Format p-values
  cex = 0.8,  # Text size
  col = "black"
)

# Make shiny object 
library(Seurat)
library(ShinyCell)
scConf = createConfig(Tcell)
makeShinyApp(seu, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell Quick Start") 
rsconnect::setAccountInfo(name='example_name',
                          token='fill-in_token_string',
                          secret='fill-in_secret_string/')


setwd('/path')
makeShinyApp(Tcell, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell-Quick-Start-Human-Tcell") 
rsconnect::setAccountInfo(name='example_name',
                          token='fill-in_token_string',
                          secret='fill-in_secret_string/')

rsconnect::deployApp('/path/shinyApp')

