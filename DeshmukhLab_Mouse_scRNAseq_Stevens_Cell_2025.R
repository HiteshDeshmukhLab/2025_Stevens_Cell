#Codebase for murine scRNAseq analysis in Stevens et al., Microbiota-derived inosine programs protective of CD8+ T cell responses against influenza in newborns, Cell, 2025

#Main scRNAseq analysis
rm(list=ls())
library(tidyverse)
library(seuratTools)
library(Seurat)
library(Matrix)
library(reshape2)
setwd('/path/codename_SCENIC')
seu <- readRDS('/path/lung_Tcells.RDS')

# 23280 features across 11104 samples within 4 assays 
head(seu[[]])
seu 

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

# seu <- NormalizeData(seu)

seu <- FilterGenes(object = seu, min.value = 1, min.cells = 100)

seu # 1282 features across 11104 samples within 2 assays 

# Calculate Cell-Cycle with Seurat, the list of genes comes with Seurat (only for human)
s.genes <- (cc.genes$s.genes)
g2m.genes <- (cc.genes$g2m.genes)
Genes <-rownames(seu)
s.genes  <- s.genes[s.genes  %in% Genes]
g2m.genes <- g2m.genes[g2m.genes %in% Genes]
seu <- CellCycleScoring(seu, s.features = s.genes,g2m.features = g2m.genes,set.ident = TRUE)
saveRDS(seu, file = 'seu.all.RDS')

# Prepare object.list for integration
DefaultAssay(seu) <- 'RNA'
unique(seu$group)
object.list <- SplitObject(seu, split.by = 'group')
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
# Integrate data sets
seu <- IntegrateData(anchorset = anchors, dims = 1:20)
head(seu[[]])

# Normal workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.1)
seu <- RunUMAP(object = seu, dims = 1:10)
Idents(seu) <- 'seurat_clusters'
DimPlot(seu, label = T, repel = T, split.by = 'group')+NoLegend()
DimPlot(seu, label = T, repel = T)+NoLegend()

# Get immune markers  
Immune.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(Immune.markers)

seu <- subset(seu, idents = c('5'), invert =T)

# Normal workflow
DefaultAssay(seu) <- 'RNA'
seu <- ScaleData(object = seu)
seu <- FindVariableFeatures(seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:5)
seu <- FindClusters(seu, resolution = 0.2)
seu <- RunUMAP(object = seu, dims = 1:10)
Idents(seu) <- 'seurat_clusters'
DimPlot(seu, label = T, repel = T, split.by = 'group')+NoLegend()
DimPlot(seu, label = T, repel = T)+NoLegend()

Immune.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
View(Immune.markers)
write_delim(Immune.markers, file = "Immune.markers.tsv", delim = '\t')
table(seu$seurat_clusters)
# 0    1    2    3    4    5 
# 3121 2176 1389 1151  771  695 
# 4 is ILC2

# Remove ILC2
seu <- subset(seu, idents = c('4'), invert =T)

DefaultAssay(seu) <- 'RNA'
seu <- ScaleData(object = seu)
seu <- FindVariableFeatures(seu)
seu <- RunPCA(object = seu)
seu <- FindNeighbors(seu, dims = 1:5)
seu <- FindClusters(seu, resolution = 0.2)
seu <- RunUMAP(object = seu, dims = 1:10)
Idents(seu) <- 'current.id'
DimPlot(seu, label = T, repel = T, split.by = 'group')+NoLegend()
DimPlot(seu, label = T, repel = T)+NoLegend()

# Collapse cluster 2 and 5 into one cluster 
cluster.ids <- c('0','1','2','3','4','2') 
names(cluster.ids) <- levels(seu)
head(seu[[]])
seu <- RenameIdents(seu, cluster.ids)
seu[["current.id"]] <- Idents(object = seu)
unique(seu$current.id)

# This is Final UMAP, save 
saveRDS(seu, file = "/path/Tcell.final.umap.RDS")
Immune.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
View(Immune.markers)

# Return top 5 markers for cluster specified 'x'
gen_marker_table <- function(x){
  Immune.markers[Immune.markers$cluster == x, ] %>%
    head(n=5)
}

# Create a data frame of results for clusters 0-5
top5_markers <- map_dfr(0:5, gen_marker_table)
View(top5_markers)
markers <- read.csv(file = 'Final_Cluster_Markers.csv')
markers <- markers[!duplicated(markers$Gene), ]

DefaultAssay(seu) <- 'RNA'
DimPlot(seu, label = T, split.by = 'group')
FeaturePlot(seu, features = CD8.naive, label = T) # 1 is Naive CD8 T cell
FeaturePlot(seu, features = c('Cd4', 'Cd40lg', 'Foxp3'), label = T) # 4 is Treg
FeaturePlot(seu, features = c('Lag3', 'Prf1', 'Cd8a', 'Havcr2', 'Gzmb', 'Nkg7', 'Cd8b1', 'Mki67'), label =T) # 0 is effector 
FeaturePlot(seu, features = c('Foxp3', 'Il2ra', 'Ctla4', 'Tigit', 'Icos', 'Batf', 'Cd4'), label =T)
FeaturePlot(seu, features = c('Nr4a1', 'Batf', 'Sox4', 'Atf', 'Relb', 'Clock', 'Stamf6', 'Ltb', 'Id3', 'Sell'), label =T)FeaturePlot(seu, features = c('Adgre1', 'Csf1r', 'Lyz2'), label =T) # 6 and 7  is macrophage cluster, needs to be eliminated
FeaturePlot(seu, features = 'Foxp3', label = T) # 4 is Treg cells
FeaturePlot(seu, features = c('Pdcd1', 'Gzmb', 'Prf1', 'Lag3'), label = T) # 3 is early-activated CD8 Tcells
FeaturePlot(seu, features = c('Tcf7', 'Gzmk', 'Il7r', 'Sell', 'Pdcd1', 'Tnfrsf9'), label = T) # 0  is CD8 Mempory Tcells
FeaturePlot(seu, features = c('Prf1', 'Ncr1', 'Il18rap'), label = T) # 3 is NK cell 

# Remove clusters 3, 6 and 7
seu <- subset(seu, idents = c('3','6','7'), invert =T)

FeaturePlot(seu, features = c('Slamf6', 'Id3', 'Pdcd1', 'Cx3cr1'), label =T)
FeaturePlot(seu, features = c('Ccr7', 'Lef1', 'Sell'), label =T)
FeaturePlot(seu, features = c('Il7r', 'Sidt1', 'Lcn4', 'Gpr183', 'Ccr7', 'Lef1'), label =T)
FeaturePlot(seu, features = c('Lgals3', 'Hist1h2ap', 'Ube2c', 'Lgals1', 'S100a4', 'S100a6',  'Mki67', 'Anxa2'), label =T) # Texint (From Daniel et al Nat Imm) # Cluster 2 
FeaturePlot(seu, features = c('Cxcr6','Cd44','Itga4', 'Itgb7', 'Itgb1'), label =T) # Lung residency, cluster 3
FeaturePlot(seu, features = c('Bcl2',  'Bcl2a1b', 'Mki67'), label =T)

seu.0 <- subset(seu, idents = c('0'))
VlnPlot(seu.0, features = c('Nfil3', 'Adora2a'), split.by = 'group', cols = alpha(default_20, 0.5))
RidgePlot(seu.0, features = c('Nfil3'), group.by = 'group', cols = alpha(default_20, 0.5))

library(scales)
library(tidyverse)
default_20 <- c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')
tol_high_contrast_palette <- c("#BB5566", "#004488")
tol_vibrant_palette <- c("#0077BB", "#33BBEE", "#009988",
                                  "#EE7733", "#CC3311", "#EE3377",
                                  "#BBBBBB")
tol_muted_palette <- c("#332288", "#ff7f0e", "#999933", "#CC6677", "#882255", "#AA4499",'#aa40fc')
DimPlot(seu, label = F,cols = alpha(tol_muted_palette, 0.5), pt.size = 1, shuffle = T, repel = T)

FeaturePlot(seu, features = c('Klrg1', 'Ly6c2', 'Klrb1c', 'Klrd1', 'S1pr1'), label =T) # Cluster 0 
FeaturePlot(seu, features = c('Tcf7', 'Pdcd1', 'Ctla4', 'Tox', 'Havcr2', 'Mki67', 'Gzmb'), label =T) 
FeaturePlot(seu, features = c('Pdcd1', 'Ctla4', 'Lag3', 'Tigit', 'Havcr2', 'Tox'), label =T) # 0 is CD8 Terminally differenriated efector cell. 
FeaturePlot(seu, features = c('Gzmb', 'Gzmk', 'Tcf7', 'Gzmm', 'Cxcr3', 'Lcn4', 'Pdcd1', 'Sidt1'), label =T) # 3 - TEM. 
FeaturePlot(seu, features = c('Slamf6', 'Cxcl10', 'Id3', 'Tcf7', 'Pdcd1'),label =T)
FeaturePlot(seu, features = c('Cx3cr1', 'Cxcr6', 'Havcr2', 'Pdcd1'),label =T)
FeaturePlot(seu, features = c('Nfil3', 'Foxo1', 'Lef1'),label =T, min.cutoff = 0.0001)

# Plot source of variations between the samples, IDs and groups 
DefaultAssay(seu) <- 'RNA'
bulk <- AggregateExpression(seu, return.seurat = TRUE, group.by = c('current.id'))
bulk <- NormalizeData(bulk)%>% FindVariableFeatures() %>% ScaleData()%>% RunPCA(npcs = 10)
bulk$group <- sapply(strsplit(Cells(bulk), split = '_'), '[',1)
bulk$celltype <- sapply(strsplit(Cells(bulk), split = '_'), '[',2)
bulk$id <- paste0(bulk$group,bulk$celltype)
A <- DimPlot(bulk, group.by = 'celltype', label = TRUE,  reduction = 'pca', pt.size = 7, cols = alpha(default_20, 0.5))
B <- DimPlot(bulk, group.by = 'group', label = TRUE,  reduction = 'pca', cols = alpha(c("#332288", "#999933"), 0.5), pt.size = 7)
C <- DimPlot(bulk, group.by = 'id', label = TRUE,  reduction = 'pca', cols = alpha(default_20, 0.5), pt.size = 7)
A+B+C # This is PCA plot 
bulk <- AggregateExpression(seu, group.by = c('current.id', 'group'))
bulk <- AggregateExpression(seu, group.by = c('current.id'))
# sig_counts <- bulk$RNA[rownames(bulk$RNA) %in% top5_markers$gene, ]
sig_counts.F <- bulk$RNA[rownames(bulk$RNA) %in% markers,]

## Set a palette
library(RColorBrewer)
library(pheatmap)
heat_colors <- rev(brewer.pal(11, "RdBu"))

## Run pheatmap using the metadata data frame for the annotation
A <- pheatmap(sig_counts.F, 
         color = heat_colors, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = T,show_colnames = T,
         border_color = NA, 
         fontsize = 10,
         scale="row",
         fontsize_row = 6, 
         height = 20)
# For complex heatmap
scaled_mat = t(scale(t(sig_counts.F)))
Heatmap(scaled_mat, col = heat_colors)

# DotPlot 
library(viridis)
library(scCustomize)
DefaultAssay(seu) <- 'RNA'
DotPlot(seu, features = markers.dp, dot.scale = 6) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  RotatedAxis() +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

# Now lets get DEG for each clusters 
Idents(seu) <- 'group'
group <- FindMarkers(object = seu, ident.1 = "NABX", ident.2 = "ABX", test.use = 'wilcox', logfc.threshold = 0.25)

bulk <- AggregateExpression(seu, group.by = c('seurat_clusters', 'group'))
sig_counts <- bulk$RNA[rownames(bulk$RNA) %in% rownames(group), ]
## Set a palette
heat_colors <- rev(brewer.pal(11, "RdBu"))
heat_colors <- viridis(n = 256, alpha = 0.7, 
                begin = 0, end = 1, option = "viridis")

## Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_counts, 
         color = heat_colors, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = T,show_colnames = T,
         scale = 'row',
         border_color = NA, 
         fontsize = 10,
         fontsize_row = 6, 
         height = 20)  

# Ridgeplots
VlnPlot(seu, features = c('Cd44', 'Foxp1', 'Il4ra', 'Lgals1', 'Lgals3'), split.by = 'group', pt.size = 0)

# Monocle 3
# Pseudotime analysis
library(SeuratWrappers)
library(monocle3)

DefaultAssay(seu) <- "RNA"
library(monocle3)
library(SeuratWrappers)
cds <- as.cell_data_set(seu)
cds <- cluster_cells(cds)
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", group_label_size = 3, show_trajectory_graph = TRUE)
plot_cells(cds, color_cells_by = "seurat_clusters", show_trajectory_graph = T)

cds <- estimate_size_factors(cds)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "current.id", label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 2, show_trajectory_graph = FALSE, cell_size = 0.8, alpha = 0.6)
rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F, label_roots = F, 
           show_trajectory_graph = F, cell_size = 0.6, alpha = 0.6)

# DEG as per monocle 3
marker_test_res <- top_markers(cds, group_cells_by="current.id", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(4, pseudo_R2)

rowData(cds)$gene_short_name <- row.names(rowData(cds))
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="current.id",
                    ordering_type="maximal_on_diag",
                    max.size=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="current.id",
                    ordering_type="cluster_row_col",
                    max.size=8
                   )

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
Idents(integrated.sub) <- "treat"
ABX <- subset(integrated.sub, ident = "ABX")
FeaturePlot(ABX, "monocle3_pseudotime", pt.size = 1.5)
NABX <- subset(integrated.sub, ident = "NABX")
FeaturePlot(NABX, "monocle3_pseudotime", pt.size = 1.5)
FeaturePlot(integrated.sub, "monocle3_pseudotime", pt.size = 1.5)

# Workflow for identifying genes which change as function of pseudotime #
cds <- estimate_size_factors(cds)
cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
head(cds_graph_test_results)
deg_ids <-rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.01))
deg_ids <-rownames(subset(cds_graph_test_results, q_value < 0.01))
view(deg_ids)
write.table(deg_ids, file = "deg_ids.csv", sep =",")
plot_cells(cds, genes = head(deg_ids), show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE)

# Collect the trajectory variable genes into modules
gene_modules <- find_gene_modules(cds[deg_ids,], resolution = c(10^seq(-2,-1)))
View(gene_modules)
write.table(gene_modules, file = "Path_to_save/Neutrophil/genes_changing_as_pseudotime.csv", sep=",") # not used in the manuscript 

#Visualize the GENES that change as function of pseudotime #
plot_cells(cds, genes=gene_modules,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE, cell_size = 1.5)

module4 <- filter(gene_modules, module ==4)
plot_cells(cds, genes=module4,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE, cell_size = 2)

module3 <- filter(gene_modules, module ==3)
plot_cells(cds, genes=module3,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE, cell_size = 2)

Interesting_genes.2 <- c('Lef1', 'Tcf7', 'Sell')
Interesting_cds <- cds[rowData(cds)$gene_short_name %in% Interesting_genes.2]

plot_genes_in_pseudotime(Interesting_cds, min_expr= 0.01, cell_size = 0.75)

#Codebase for murine NFIL3 ChIPseq analysis in Stevens et al., Microbiota-derived inosine programs protective of CD8+ T cell responses against influenza in newborns (working title), Cell, 2025
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DiffBind)
library(DESeq2)
library(tidyverse)
library(GenomicRanges)
library(rjson)
library(tidyverse)
setwd("/path")

samples <- read.delim('NFIL3_1.csv', sep = ",")
names(samples)
samples
NFIL3 <- dba(sampleSheet=samples)
NFIL3
plot(NFIL3)

# Calculate affinity matrix 
NFIL3 <- dba.count(NFIL3)
NFIL3

# Counting reads
info <- dba.show(NFIL3)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
plot(NFIL3)
dev.off()
# Normalizing the data
NFIL3 <- dba.normalize(NFIL3)
norm <- dba.normalize(NFIL3, bRetrieve=TRUE)
normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors, NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs
# Establishing a model design and contrast
NFIL3 <- dba.contrast(NFIL3, minMembers=2, reorderMeta=list(Condition="NoABX"))
# NoABX condition is  the baseline, so it will be in the denominator of the default contrast.
NFIL3

# Performing the differential analysis
NFIL3 <- dba.analyze(NFIL3)
dev.off()
par(mar = c(1, 1, 1, 1))
plot(NFIL3, contrast=1)

# Retrieving the differentially bound sites
NFIL3.DB <- dba.report(NFIL3)
NFIL3.DB 

# Plotting in DiffBind
dba.plotVolcano(NFIL3, dotSize = 0.5) 
# Venn Plot
dba.plotVenn(NFIL3, contrast=1, bDB=TRUE,bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dba.plotPCA(NFIL3,contrast=1,label=DBA_ID)
# MA plot 
dba.plotMA(NFIL3)

# Heatmaps
corvals <- dba.plotHeatmap(NFIL3)
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
dev.off()
readscores <- dba.plotHeatmap(NFIL3, contrast=1, correlations=FALSE, scale="row", colScheme = hmap)

# Profiling and Profile Heatmaps
BiocManager::install("profileplyr")
profiles <- dba.plotProfile(NFIL3)
dba.plotProfile(profiles)

# Merging all samples in a contrast condition
merge=c(DBA_TISSUE, DBA_REPLICATE)
profiles <- dba.plotProfile(NFIL3,merge=c(DBA_TISSUE, DBA_REPLICATE))
dba.plotProfile(profiles)

dba.plotHeatmap(NFIL3)

# Retrieve affinity matrix 
NFIL3.DB <- dba.report(NFIL3)
NFIL3.DB
write.csv(NFIL3.DB, file ="DE.Peaks.csv")

# Create bed files for each keeping only significant peaks (p < 0.05)
out <- as.data.frame(NFIL3.DB)

NABX_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)
# Write to file
write.table(NABX_enrich, file="NABX_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

ABX_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)
# Write to file
write.table(ABX_enrich, file="ABX_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Annotate Peaks with homer 
# changeNewLine.pl DE_Peaks_bed.txt
# annotatePeaks.pl  DE_Peaks_bed.txt mm10 > DE_anno.txt


# For visualizing the tracks
# bamCoverage -b bowtie2_on_data_9.bam \
# -o visualization/bigWig/bowtie2_on_data_9.bw \
# --binSize 20 \
# --normalizeUsing BPM \
# --smoothLength 60 \
# --extendReads 150 \
# --centerReads \
# -p 6 2> ../logs/Nanog_rep2_bamCoverage.log

# ChipAnno 
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
samplefiles <- list.files("/path/peaks", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
print(samplefiles)
peak <- readPeakFile(samplefiles[[4]])
peak
# ChIP peaks coverage plot
covplot(peak, weightCol="V5")
# Saved as PDF 

# Profile of ChIP peaks binding to TSS regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
head(promoter)
tagMatrix <- getTagMatrix(peak, windows=promoter)
head(tagMatrix)
# Heatmap of ChIP binding to TSS regions
tagHeatmap(tagMatrix)

peakHeatmap(samplefiles[[4]], TxDb=txdb, upstream=3000, downstream=3000)
# Saved heatmap as PDF

# Peak Annotation
peakAnno <- annotatePeak(samplefiles[[4]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
# Visualize Genomic Annotation
plotAnnoPie(peakAnno)
head(as.data.frame(peakAnno))

edb <- EnsDb.Mmusculus.v79
seqlevelsStyle(edb) <- "UCSC"
peakAnno.edb <- annotatePeak(samplefiles[[4]], tssRegion=c(-3000, 3000),
                             TxDb=edb, annoDb="org.Mm.eg.db")
head(as.data.frame(peakAnno.edb))

# Functional enrichment analysis
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId, organism = "mouse")
head(pathway1, 2)
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)

# Additional code used for NFIL3 pseudotime analysis
Idents(seu) <- 'seurat_clusters'
DotPlot(seu, features = 'Nfil3')

DimPlot(seu)
cds <- as.cell_data_set(seu)
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

# Marker
rowData(cds)$gene_short_name <- row.names(rowData(cds))
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.2) %>%
  group_by(cell_group) %>%
  top_n(4, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- seu@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seu@reductions$umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, show_trajectory_graph = FALSE,
                                 group_label_size = 0, cell_size = 0.75,
                                 alpha = 0.5) + theme(legend.position = "right")
cluster.before.traj
cds <- learn_graph(cds, use_partition = T)

# Ordering cells 
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           label_roots = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 0.75,
           alpha = 0.5)

cds <- estimate_size_factors(cds)
cds_graph_test_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
cds_graph_test_results %>% arrange(q_value) %>% filter(status == "OK") %>% head()

# Cell ordered by Monocle3
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()

#Add pseudotime values into the seuratobject
seu$pseudotime <- pseudotime(cds)
FeaturePlot(seu, features = "pseudotime")

RidgePlot(seu, features = c('Lef1', 'Cd44', 'Havcr2', 'Pdcd1'), sort = T, idents = c("1", "4", "3", "0", "5", "2"))
RidgePlot(seu, features = c('Lef1', 'Tcf7', 'Havcr2', 'Pdcd1'), sort = T, idents = c("1", "0"))

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Lef1", "Tcf7", "Zeb2", "Slamf6", "Nfil3", "Sell"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


head(cds_graph_test_results)
deg_ids <-rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.1))
deg_ids <-rownames(subset(cds_graph_test_results, q_value < 0.1))
view(deg_ids)
write.table(deg_ids, file = "deg_ids.csv", sep =",")
plot_cells(cds, genes = c('Pdcd1', 'Cd44'), show_trajectory_graph = TRUE, label_cell_groups = FALSE, label_leaves = FALSE)

# End code
