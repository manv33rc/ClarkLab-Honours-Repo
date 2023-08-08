## Integrating Day 25 and Day 55 datasets with Harmony

# Load libraries
library(harmony)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)
library(reticulate)

# Set seed for reproducibility
set.seed(4242)

# Import our feature_ID_converter python script
source_python("/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/feature_ID_converter.py")
gtfFilePath <- "/data/gpfs/projects/punim0646/manveer/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

# Load cell cyle markers obtained from : https://hbctraining.github.io/scRNA-seq_online/lessons/06_SC_SCT_normalization.html
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")

# Create empty lists to store cluster resolution figures for integrated short read and long read Seurat objects
SR.cluster_resolution.figs <- list()
LR.cluster_resolution.figs <- list()

regressCellPhase <- TRUE

## Read in filtered seurat objects for both samples with short reads and long reads---------------
sr.day25.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/SR_DAY25-QCed.rds"
sr.day55.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/SR_DAY55-QCed.rds"

lr.day25.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/LR_DAY25-QCed.rds"
lr.day55.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/LR_DAY55-QCed.rds"

lr.day25.isoLevel.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/LR_DAY25-isoLvl-QCed.rds"
lr.day55.isoLevel.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/LR_DAY55-isoLvl-QCed.rds"

# Use filepaths to .RDS objects, so that we have the seurat objects we need to work with in this environment
short.read.day25 <- readRDS(file = sr.day25.filepath)
short.read.day55 <- readRDS(file = sr.day55.filepath)
long.read.day25 <- readRDS(file = lr.day25.filepath)
long.read.day55 <- readRDS(file = lr.day55.filepath)

long.read.day25.IsoLvl <- readRDS(file = lr.day25.isoLevel.filepath)
long.read.day55.IsoLvl <- readRDS(file = lr.day55.isoLevel.filepath)

# Merge the day 25 and 55 seurat objects for short reads and long reads into individual seurat objects------------------
short.reads.combined <- merge(x = short.read.day25,
                              y = short.read.day55,
                              project = "Combined.SR.Samples")
long.reads.combined <- merge(x = long.read.day25,
                              y = long.read.day55,
                              project = "Combined.LR.Samples")
long.reads.IsoLvl.combined <- merge(x = long.read.day25.IsoLvl,
                                    y = long.read.day55.IsoLvl,
                                    project = "Combined.LR.Samples.IsoLvl")


### Cell number, feature number, UMI count comparisons between samples----------------------
# Code sourced from : https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

## Compare Cell numbers
compareCellNumbers <- function(seurat.obj, title = ""){
  plt <- seurat.obj@meta.data %>%
    ggplot(aes(x = seurat.obj$orig.ident, fill = seurat.obj$orig.ident)) +
    geom_bar() + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle(paste0("Number of Cells (",title,")")) +
    xlab("Sample")
  
  plt
}
cell.num.bar.plt.SR <- compareCellNumbers(short.reads.combined, title = "Integrated Short Reads")
cell.num.bar.plt.LR <- compareCellNumbers(long.reads.combined, title = "Integrated Long Reads")

## Compare UMI counts (transcripts) per cell
compareReadCounts <- function(seurat.obj, title = ""){
  plt <- seurat.obj@meta.data %>% 
    ggplot(aes(color = seurat.obj$orig.ident, x = seurat.obj$nCount_RNA,
               fill = seurat.obj$orig.ident)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell Density") +
    xlab("nUMI") +
    geom_vline(xintercept = 1000) +
    labs(title = paste0("UMI counts per cell (",title,")"))
  
  plt
}
UMI.density.plt.SR <- compareReadCounts(short.reads.combined, "Short Reads")
UMI.density.plt.LR <- compareReadCounts(long.reads.combined, "Long Reads")

## Visualize the distribution of genes detected per cell to compare short reads and long reads
compareFeaturesDetected <- function(seurat.obj, title = "", featureType = ""){
  seurat.obj@meta.data %>% 
    ggplot(aes(color = seurat.obj$orig.ident, x = seurat.obj$nFeature_RNA,
               fill = seurat.obj$orig.ident)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = 1000) +
    ylab("Cell Density") +
    xlab(paste0("Number of ", featureType)) +
    labs(title = paste0(featureType, " detected per cell (",title,")"))
}
Gene.density.plt.SR <- compareFeaturesDetected(short.reads.combined, "Short Reads", "Genes")
Gene.density.plt.LR <- compareFeaturesDetected(long.reads.combined, "Long Reads", "Genes")
Iso.density.plt.LR <- compareFeaturesDetected(long.reads.IsoLvl.combined, "Long Reads", "Isoforms")

## Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
compareUMIGeneCorrelations <- function(seurat.obj, title = "", xint = 1000, yint = 1000){
  seurat.obj@meta.data %>% 
    ggplot(aes(x=seurat.obj$nCount_RNA, y=seurat.obj$nFeature_RNA, color=seurat.obj$percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = xint) +
    geom_hline(yintercept = yint) +
    xlab("nUMI") +
    ylab("nGenes") +
    labs(title = paste0("UMI counts and Genes per Cell Correlation\n(", title,")"))
}
LR.UMI.Gene.Association.plt.BEFORE.FILTER <- compareUMIGeneCorrelations(long.reads.combined,
                                                                        "Integrated Long Reads",
                                                                        xint = 7400, yint = 2600)
SR.UMI.Gene.Association.plt.BEFORE.FILTER <- compareUMIGeneCorrelations(short.reads.combined,
                                                                        "Integrated Short Reads",
                                                                        xint = 10500, yint = 4800)
LR.UMI.Gene.Association.plt.BEFORE.FILTER | SR.UMI.Gene.Association.plt.BEFORE.FILTER

## Filter out low quality reads using selected thresholds
# After looking at previous Association Plots, Both integrated Seurat objects require extra filtering
# Poor quality cells likely to have low genes and UMIs per cell (bottom left plot quadrant). 
# Good cells will generally exhibit both higher number of genes AND UMIs per cell
short.reads.combined.filtered <- subset(x = short.reads.combined, 
                          subset = (nCount_RNA >= 10500) &
                            (nFeature_RNA >= 4800))
long.reads.combined.filtered <- subset(x = long.reads.combined,
                              subset = (nCount_RNA >= 7400 &
                                        nFeature_RNA >= 2600))
LR.UMI.Gene.Association.plt.AFTER.FILTER <- compareUMIGeneCorrelations(long.reads.combined.filtered,
                                                                       "Filtered Integrated Long Reads",
                                                                       xint = 7000,
                                                                       yint = 2500)
SR.UMI.Gene.Association.plt.AFTER.FILTER <- compareUMIGeneCorrelations(short.reads.combined.filtered,
                                                                       "Filtered Integrated Short Reads",
                                                                       xint = 10000,
                                                                       yint = 5400)

# Evaluate the effects of filtering with our different plots----------------------------------------
LR.UMI.Gene.Association.plt.AFTER.FILTER | SR.UMI.Gene.Association.plt.AFTER.FILTER
LR.UMI.Gene.Association.plt.BEFORE.FILTER | LR.UMI.Gene.Association.plt.AFTER.FILTER
SR.UMI.Gene.Association.plt.BEFORE.FILTER | SR.UMI.Gene.Association.plt.AFTER.FILTER

UMI.density.plt.SR.AfterFilter <- compareReadCounts(short.reads.combined.filtered, "Filtered Short Reads")
UMI.density.plt.SR | UMI.density.plt.SR.AfterFilter
Gene.density.plt.SR.AfterFilter <- compareFeaturesDetected(short.reads.combined.filtered, "Filtered Short Reads", "Genes")
Gene.density.plt.SR | Gene.density.plt.SR.AfterFilter

UMI.density.plt.LR.AfterFilter <- compareReadCounts(long.reads.combined.filtered, "Filtered Long Reads")
UMI.density.plt.LR | UMI.density.plt.LR.AfterFilter
Gene.density.plt.LR.AfterFilter <- compareFeaturesDetected(long.reads.combined.filtered, "Filtered Long Reads", "Genes")
Gene.density.plt.LR | Gene.density.plt.LR.AfterFilter

cell.num.bar.plt.SR.FILTERED <- compareCellNumbers(short.reads.combined.filtered, title = "Filtered Short Reads")
cell.num.bar.plt.LR.FILTERED <- compareCellNumbers(long.reads.combined.filtered, title = "Filtered Long Reads")
cell.num.bar.plt.LR | cell.num.bar.plt.LR.FILTERED
cell.num.bar.plt.SR | cell.num.bar.plt.SR.FILTERED

## Make sure that the cells in the filtered gene-level Long Read Seurat object match those found in the Isoform-level object
long.reads.IsoLvl.combined <- subset(long.reads.IsoLvl.combined, 
                                     cells = row.names(long.reads.combined.filtered@meta.data))
#identical(row.names(long.reads.IsoLvl.combined@meta.data),row.names(long.reads.combined.filtered@meta.data))



### Assign a score to each cell based on its expression of G2/M and S phase markers for long reads and short reads--------------
### These scores will be used later with PCA to determine whether cell cycle is a major source of variation
short.reads.combined.filtered <- CellCycleScoring(short.reads.combined.filtered,
                                         g2m.features = g2m_genes,
                                         s.features = s_genes)

# convert our cell cycle markers to their respective IDs to use with our Long-read gene lvl seurat object
g2m_gene.IDs <- ListofGeneNamestoIDs(g2m_genes, gtfFilePath)
s_gene.IDs <- ListofGeneNamestoIDs(s_genes, gtfFilePath)
long.reads.combined.filtered <- CellCycleScoring(long.reads.combined.filtered,
                                        g2m.features = g2m_gene.IDs,
                                        s.features = s_gene.IDs)


### perform standard seurat workflow steps for short read data------------------------------------
str(short.reads.combined.filtered) # They haven't been scaled, and variable features slot is empty
short.reads.combined.filtered[['RNA']]@data@x # Our counts have already been normalized

short.reads.combined.filtered <- FindVariableFeatures(short.reads.combined.filtered,
                                             selection.method = 'vst',
                                             nfeatures = 2000)

if(regressCellPhase == TRUE){
  short.reads.combined.filtered <- ScaleData(short.reads.combined.filtered, vars.to.regress = c("S.Score", "G2M.Score"))
} else{
  short.reads.combined.filtered <- ScaleData(short.reads.combined.filtered)
}

short.reads.combined.filtered <- RunPCA(short.reads.combined.filtered,
                               features = VariableFeatures(object = short.reads.combined.filtered))
elbow.plt.combined.short.reads <- ElbowPlot(short.reads.combined.filtered) +
  labs(title = "Elbow Plot - Integrated Short Reads")
elbow.plt.combined.short.reads
short.reads.combined.filtered <- RunUMAP(short.reads.combined.filtered, dims = 1:12, reduction = 'pca') # Choosing 12 PCs based on elbow plot

before.harmony.shortReads <- DimPlot(short.reads.combined.filtered, reduction = 'umap', group.by = 'orig.ident')
before.harmony.shortReads

### perform standard seurat workflow steps for long read data (Gene Lvl)----------------------------------------
long.reads.combined.filtered <- FindVariableFeatures(long.reads.combined.filtered,
                                             selection.method = 'vst',
                                             nfeatures = 2000)
if(regressCellPhase == TRUE){
  long.reads.combined.filtered <- ScaleData(long.reads.combined.filtered, vars.to.regress = c("S.Score", "G2M.Score"))
} else{
  long.reads.combined.filtered <- ScaleData(long.reads.combined.filtered)
}
long.reads.combined.filtered <- RunPCA(long.reads.combined.filtered,
                               features = VariableFeatures(object = long.reads.combined.filtered))
elbow.plt.combined.long.reads <- ElbowPlot(long.reads.combined.filtered) + 
  labs(title = "Elbow Plot - Integrated Long Reads\n(Gene Level)")
elbow.plt.combined.long.reads
long.reads.combined.filtered <- RunUMAP(long.reads.combined.filtered, dims = 1:11, reduction = 'pca') # Choosing 11 PCs based on elbow plot

before.harmony.longReads <- DimPlot(long.reads.combined.filtered, reduction = 'umap', group.by = 'orig.ident')
before.harmony.longReads



### perform standard seurat workflow steps for long read data (Isoform Lvl)------------------------------------
long.reads.IsoLvl.combined <- FindVariableFeatures(long.reads.IsoLvl.combined,
                                             selection.method = 'vst',
                                             nfeatures = 2000)
long.reads.IsoLvl.combined <- ScaleData(long.reads.IsoLvl.combined)
long.reads.IsoLvl.combined <- RunPCA(long.reads.IsoLvl.combined,
                               features = VariableFeatures(object = long.reads.IsoLvl.combined))
elbow.plt.combined.long.reads.IsoLvl <- ElbowPlot(long.reads.IsoLvl.combined) +
  labs(title = "Elbow Plot - Integrated Long Reads\n(Isoform Level)")
elbow.plt.combined.long.reads.IsoLvl
long.reads.IsoLvl.combined <- RunUMAP(long.reads.IsoLvl.combined, dims = 1:12, reduction = 'pca') # Choosing 12 PCs based on elbow plot

before.harmony.longReads.IsoLvl <- DimPlot(long.reads.IsoLvl.combined, reduction = 'umap', group.by = 'orig.ident')
before.harmony.longReads.IsoLvl








### Run Harmony to integrate our short read data across the two maturation timepoints------------------
# This will overlay similar cells across samples together
short.reads.harmony <- short.reads.combined.filtered %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = TRUE)
short.reads.harmony@reductions

SR.harmony.embed <- Embeddings(short.reads.harmony, reduction = 'harmony')
SR.harmony.embed[1:20][1:20]

# Do UMAP and clustering using **Harmony embedding values instead of PCA**
short.reads.harmony <- short.reads.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:12) %>%  # Using the same number of PCs as with RunPCA() earlier
  FindNeighbors(reduction = "harmony", dims = 1:12) %>% 
  #FindClusters(resolution = c(0.3, 0.5, 0.7, 0.9, 1.1))
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))

# SET SHORT READ CLUSTER IDENTITY HERE - BASED ON CLUSTER RESOLUTION FIGS
short.reads.harmony$seurat_clusters <- short.reads.harmony$"RNA_snn_res.0.3"
#View(short.reads.harmony@meta.data)

UMAP.integrated.short.reads <- DimPlot(short.reads.harmony, reduction = 'umap', group.by = 'orig.ident')

# Save cluster resolution plots to determine the appropriate cluster resolution
SR.cluster_resolution.figs[[1]] <- DimPlot(short.reads.harmony, group.by = "RNA_snn_res.0.1", label = TRUE)
SR.cluster_resolution.figs[[2]] <- DimPlot(short.reads.harmony, group.by = "RNA_snn_res.0.2", label = TRUE)
SR.cluster_resolution.figs[[3]] <- DimPlot(short.reads.harmony, group.by = "RNA_snn_res.0.3", label = TRUE)
SR.cluster_resolution.figs[[4]] <- DimPlot(short.reads.harmony, group.by = "RNA_snn_res.0.4", label = TRUE)
SR.cluster_resolution.figs[[5]] <- DimPlot(short.reads.harmony, group.by = "RNA_snn_res.0.5", label = TRUE)

UMAP.integrated.short.reads <- DimPlot(short.reads.harmony, reduction = 'umap', group.by = 'orig.ident')



### Run Harmony to integrate our long read data (GENE level) across the two maturation timepoints----------------------
long.reads.harmony <- long.reads.combined.filtered %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = TRUE)
long.reads.harmony@reductions

LR.harmony.embed <- Embeddings(long.reads.harmony, reduction = 'harmony')
LR.harmony.embed[1:20][1:20]

# Do UMAP and clustering using **Harmony embedding values instead of PCA**
long.reads.harmony <- long.reads.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:12) %>%  # Using the same number of PCs as with RunPCA() earlier
  FindNeighbors(reduction = "harmony", dims = 1:12) %>% 
  FindClusters(resolution = c(0.3, 0.5, 0.7, 0.9, 1.1))

# SET LONG READ CLUSTER IDENTITY HERE - BASED ON CLUSTER RESOLUTION FIGS
long.reads.harmony$seurat_clusters <- long.reads.harmony$"RNA_snn_res.0.5"
#View(long.reads.harmony@meta.data)

UMAP.integrated.long.reads <- DimPlot(long.reads.harmony, reduction = 'umap', group.by = 'orig.ident')

# Save cluster resolution plots to determine the appropriate cluster resolution
LR.cluster_resolution.figs[[1]] <- DimPlot(long.reads.harmony, group.by = "RNA_snn_res.0.3", label = TRUE)
LR.cluster_resolution.figs[[2]] <- DimPlot(long.reads.harmony, group.by = "RNA_snn_res.0.5", label = TRUE)
LR.cluster_resolution.figs[[3]] <- DimPlot(long.reads.harmony, group.by = "RNA_snn_res.0.7", label = TRUE)
LR.cluster_resolution.figs[[4]] <- DimPlot(long.reads.harmony, group.by = "RNA_snn_res.0.9", label = TRUE)
LR.cluster_resolution.figs[[5]] <- DimPlot(long.reads.harmony, group.by = "RNA_snn_res.1.1", label = TRUE)


### Run Harmony to integrate our long read data (ISOFORM level) across the two maturation timepoints----------------------
long.reads.harmony.IsoLvl <- long.reads.IsoLvl.combined %>%
  RunHarmony(group.by.vars = 'orig.ident', plot_convergence = TRUE)
long.reads.harmony.IsoLvl@reductions

LR.harmony.embed <- Embeddings(long.reads.harmony.IsoLvl, reduction = 'harmony')
LR.harmony.embed[1:20][1:20]

# Do UMAP and clustering using **Harmony embedding values instead of PCA**
long.reads.harmony.IsoLvl <- long.reads.harmony.IsoLvl %>%
  RunUMAP(reduction = 'harmony', dims = 1:12) %>%  # Using the same number of PCs as with RunPCA() earlier
  FindNeighbors(reduction = "harmony", dims = 1:12) %>% 
  FindClusters(resolution = c(0.3, 0.5, 0.7, 0.9, 1.1))

## COPY OVER CLUSTER IDENTITIES FROM GENE LEVEL LONG READ COUNTS
# First, save UMAP with isoform level derived clusters
clusters.long.reads.IsoLvl.beforeTransfer <- DimPlot(long.reads.harmony.IsoLvl,
                                                     reduction = 'umap',
                                                     group.by = "RNA_snn_res.0.5", label = TRUE,
                                                     label.size = 5, repel = TRUE) +
  labs(title = "Long-Read Isoform-Level UMAP\nBEFORE Gene-Level Label Transfer") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))

# Transfer labels from gene level counts
long.reads.harmony.IsoLvl$seurat_clusters <- long.reads.harmony$seurat_clusters

# Save cluster and sample-labeled UMAPs
UMAP.integrated.long.reads.IsoLvl <- DimPlot(long.reads.harmony.IsoLvl, reduction = 'umap', group.by = 'orig.ident')
clusters.long.reads.IsoLvl.afterTransfer <- DimPlot(long.reads.harmony.IsoLvl, reduction = 'umap', 
                                                    group.by = 'seurat_clusters', label = TRUE,
                                                    label.size = 5, repel = TRUE) +
  labs(title = "Long-Read Isoform-Level UMAP\nAFTER Gene-Level Label Transfer") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))





### Compile additional figures showing the effects of cell cycle phase on clustering-------------------------
pc1.pc2.sr.cellphase.plt <- DimPlot(short.reads.harmony, dims = c(1,2), reduction = 'harmony', group.by = 'Phase') +
  labs(title = "Cell Phase (Short Reads)")
pc1.pc2.lr.cellphase.plt <- DimPlot(long.reads.harmony, dims = c(1,2), reduction = 'harmony', group.by = 'Phase') +
  labs(title = "Cell Phase (Long Reads)")
pc1.pc2.sr.cellphase.plt | pc1.pc2.lr.cellphase.plt

sr.CellPhase.UMAP <- DimPlot(short.reads.harmony, reduction = 'umap',
        group.by = 'Phase', label = TRUE,
        label.size = 5, repel = TRUE) +
  ggmin::theme_powerpoint() +
  labs(title = "Cell Cycle Status:\nIntegrated Short Reads") +
  theme(plot.title = element_text(size = 16))
lr.CellPhase.UMAP <- DimPlot(long.reads.harmony, reduction = 'umap',
        group.by = 'Phase', label = TRUE,
        label.size = 5, repel = TRUE) +
  ggmin::theme_powerpoint() +
  labs(title = "Cell Cycle Status:\nIntegrated Long Reads (Gene Level)") +
  theme(plot.title = element_text(size = 16))




#### EXPLORATION OF QC METRICS----------------------------------------
# UMAP of cells in each cluster by sample
DimPlot(short.reads.harmony, 
        label = TRUE, 
        split.by = "orig.ident")  + NoLegend()
View(short.reads.harmony@meta.data)
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(short.reads.harmony, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# UMAP of cells in each cluster by sample
DimPlot(long.reads.harmony, 
        label = TRUE, 
        split.by = "orig.ident")  + NoLegend()
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
FeaturePlot(long.reads.harmony, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


### Export our integrated seurat objects to use in other scripts------------------------
# Ensure chosen cluster identities are saved
Idents(long.reads.harmony) <- long.reads.harmony$seurat_clusters
Idents(long.reads.harmony.IsoLvl) <- long.reads.harmony.IsoLvl$seurat_clusters
Idents(short.reads.harmony) <- short.reads.harmony$seurat_clusters

saveRDS(long.reads.harmony, file = "INTEGRATED_LR_days25+55.rds")
saveRDS(long.reads.harmony.IsoLvl, file = "INTEGRATED_LR_days25+55_ISOLVL.rds")
saveRDS(short.reads.harmony, file = "INTEGRATED_SR_days25_55.rds")


### Inspect before and after integration UMAPS, and the cluster labeling for short reads and long reads-------------------
## Long-read Gene-level Figures
LRclusterUMAP <- DimPlot(long.reads.harmony, reduction = 'umap', 
                         group.by = 'seurat_clusters', label = TRUE,
                         label.size = 5, repel = TRUE) +
  labs(title = "Integrated Day 25 + 55\nLong-Read Gene Level Clusters\n(NPCs = 12, Resolution = 0.5)") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
LRclusterUMAP

# Before and after integration UMAP plots
before.harmony.longReads <- before.harmony.longReads +
  labs(title = "Long-Read Gene-Level UMAP\nBefore Sample Integration") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
UMAP.integrated.long.reads <- UMAP.integrated.long.reads + 
  labs(title = "Long-Read Gene-Level UMAP\nAfter Sample Integration") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
before.harmony.longReads | UMAP.integrated.long.reads

## Long-read Isoform-level Figures
clusters.long.reads.IsoLvl.beforeTransfer | clusters.long.reads.IsoLvl.afterTransfer

before.harmony.longReads.IsoLvl <- before.harmony.longReads.IsoLvl +
  labs(title = "Long-Read Isoform-Level UMAP\nBefore Sample Integration") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
UMAP.integrated.long.reads.IsoLvl <- UMAP.integrated.long.reads.IsoLvl +
  labs(title = "Long-Read Isoform-Level UMAP\nAfter Sample Integration") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
before.harmony.longReads.IsoLvl | UMAP.integrated.long.reads.IsoLvl

## Short-read Gene-level Figures
SRclusterUMAP <- DimPlot(short.reads.harmony, reduction = 'umap', 
                         group.by = 'seurat_clusters', label = TRUE,
                         label.size = 5) +
  labs(title = "Integrated Day 25 + 55\nShort-Read Clusters\n(NPCs = 12, Resolution = 0.3)") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
SRclusterUMAP

before.harmony.shortReads <- before.harmony.shortReads +
  labs(title = "Short-Read Isoform-Level UMAP\nBefore Sample Integration") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
UMAP.integrated.short.reads <- UMAP.integrated.short.reads +
  labs(title = "Short-Read Isoform-Level UMAP\nAfter Sample Integration") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
before.harmony.shortReads | UMAP.integrated.short.reads


### Save our UMAPs into a PDFs---------------------------------------------------
pdf("Harmony_Integration_UMAPs.pdf", width = 15, height = 15)

integration.compare.LRgenelvl <- grid.arrange(before.harmony.longReads, UMAP.integrated.long.reads,
             nrow = 2, ncol = 1,
             top = textGrob("Harmony Integration (Days 25 + 55)\nLong Read Gene Level Counts\n",
                          gp = gpar(fontsize = 20)))
grid.draw(integration.compare.LRgenelvl)
grid.draw(LRclusterUMAP)

integration.compare.LRisoLvl <- grid.arrange(before.harmony.longReads.IsoLvl, UMAP.integrated.long.reads.IsoLvl,
                                             nrow = 2, ncol = 1,
                                             top = textGrob("Harmony Integration (Days 25 + 55)\nLong Read Isoform Level Counts\n",
                                                            gp = gpar(fontsize = 20)))
grid.draw(integration.compare.LRisoLvl)
labeltransfer.compare.LRisolvl <- grid.arrange(clusters.long.reads.IsoLvl.beforeTransfer,
                                               clusters.long.reads.IsoLvl.afterTransfer,
                                               nrow = 2, ncol = 1,
                                               top = textGrob("Before and After Transferring Labels\nderived from Gene-Level counts\n",
                                                              gp = gpar(fontsize = 20)))
grid.draw(labeltransfer.compare.LRisolvl)

integration.compare.SR <- grid.arrange(before.harmony.shortReads, UMAP.integrated.short.reads,
                                       nrow = 2, ncol = 1,
                                       top = textGrob("Harmony Integration (Days 25 + 55)\nShort Read Counts",
                                                      gp = gpar(fontsize = 20)))
grid.draw(integration.compare.SR)
grid.draw(SRclusterUMAP)

sr.lr.cluster.compare <- grid.arrange(LRclusterUMAP, SRclusterUMAP,
                                      nrow = 2, ncol = 1,
                                      top = textGrob("Unsupervised Clustering Comparison\nbetween SR and LR (gene level)\n",
                                                     gp = gpar(fontsize = 20)))

elbow.plots <- grid.arrange(elbow.plt.combined.long.reads, elbow.plt.combined.long.reads.IsoLvl,
                            elbow.plt.combined.short.reads,
                            nrow = 2, ncol = 2)
# Close the PDF device
dev.off()



### Save our metadata comparison plots between short reads and long reads into a PDF---------------------
pdf("Metadata-SR.LR.comparison.pdf", width = 10, height = 8)

cell.num.layout <- grid.arrange(cell.num.bar.plt.SR, cell.num.bar.plt.LR,
                                nrow = 1, ncol = 2,
                                top = textGrob("Comparing Cell Numbers per Sample (after QC)\nLong Reads and Short Reads\n"))
grid.draw(cell.num.layout)

UMI.density.layout <- grid.arrange(UMI.density.plt.LR.iso, UMI.density.plt.SR,
                                   nrow = 1, ncol = 2,
                                   top = textGrob("Comparing Reads per Cell across Samples (after QC)\n"))
grid.draw(UMI.density.layout)

Feature.density.layout <- grid.arrange(Gene.density.plt.LR, Iso.density.plt.LR,
                                       Gene.density.plt.SR,
                                       nrow = 2, ncol = 2,
                                       top = textGrob("Comparing Features Detected per Cell across Samples\n(after QC)\n"))
grid.draw(Feature.density.layout)
# Close the PDF device
dev.off()


### Save our cell cycle figures for short reads and long reads into a PDF-------------------------------
Cell.cycle.layout <- grid.arrange(pc1.pc2.lr.cellphase.plt, pc1.pc2.sr.cellphase.plt,
                                  lr.CellPhase.UMAP, sr.CellPhase.UMAP,
                                  nrow = 2, ncol = 2,
                                  top = textGrob("Effects of Cell Cycle Phase on Unsupervised Clustering\n",
                                                 gp = gpar(fontsize = 15)))
ggsave("Cell-cycle-visualizations-BEFORE-REGRESSION.pdf", Cell.cycle.layout, width = 15, height = 15)


### Save all the cluster resolution figures for integrated long reads and short reads into PDFs---------------------------------------
## Generate Short read cluster resolution pdf
SR.cluster.res.layout <- grid.arrange(SR.cluster_resolution.figs[[1]], SR.cluster_resolution.figs[[2]],
                                   SR.cluster_resolution.figs[[3]], SR.cluster_resolution.figs[[4]], 
                                   SR.cluster_resolution.figs[[5]], SRclusterUMAP,
                                   nrow = 2, ncol = 3,
                                   top=textGrob("Harmony Integration - SHORT READ Cluster Resolution Figures",
                                                gp = gpar(fontsize = 15)))

ggsave("SR-cluster-res-figures-HARMONY.pdf", SR.cluster.res.layout, width = 20, height = 12)

## Generate Long read cluster resolution pdf
LR.cluster.res.layout <- grid.arrange(LR.cluster_resolution.figs[[1]], LR.cluster_resolution.figs[[2]],
                                      LR.cluster_resolution.figs[[3]], LR.cluster_resolution.figs[[4]], 
                                      LR.cluster_resolution.figs[[5]], LRclusterUMAP,
                                      nrow = 2, ncol = 3,
                                      top=textGrob("Harmony Integration - LONG READ Cluster Resolution Figures",
                                                   gp = gpar(fontsize = 15)))

ggsave("LR-cluster-res-figures-HARMONY.pdf", LR.cluster.res.layout, width = 20, height = 12)



### TOP VARIABLE GENE INFORMATION -     INTEGRATED LONG READS-------------------------------------------
# Extract the top 10 most variable features across cells
top10.variable.genes.lr <- head(VariableFeatures(long.reads.harmony), 10)
# Convert gene IDs in top10 to their respective symbols using FeatureIDtoName() function
converted_top10.lr <- sapply(top10.variable.genes.lr, function(feature_id) {
  FeatureIDtoName(feature_id, ConvType = 'g.to.gName', gtfFilePath = gtfFilePath)
})

## Plot the 2000 most variable features and label the top10 most variable transcripts
variable.feature.plot.lr <- VariableFeaturePlot(long.reads.harmony) +
  labs(title = "Top 2000 Variable Genes across Cells\n(Integrated LR)")
variable.feature.plot.lr <- LabelPoints(plot = variable.feature.plot.lr, points = top10.variable.genes.lr, labels = converted_top10.lr,
                                     repel = TRUE, xnudge = 0, ynudge = 0)
variable.feature.plot.lr

## Create feature plots showing the expression of these genes
FoI_figures <- list()
converted_top10.lr
for (feature in top10.variable.genes.lr) {
  feature_name <- converted_top10.lr[feature]
  FoI_figures[[feature]] <- FeaturePlot(long.reads.harmony, reduction = 'umap',
                                        features = feature,
                                        order = TRUE) + ggmin::theme_powerpoint() + labs(title = feature_name)
}

grid.layout.FoI <- grid.arrange(grobs = FoI_figures, nrow = 2, ncol = 5,
                                top=textGrob("Expression Levels of Top 10 Most Variable Genes\n(Integrated Long Reads)\n",
                                             gp = gpar(fontsize = 12, fontface = "bold")))
ggsave("Integrated-LR-top10-genes.pdf", grid.layout.FoI, width = 22, height = 10)

### TOP VARIABLE GENE INFORMATION -     INTEGRATED SHORT READS-------------------------------------------
# Extract the top 10 most variable features across cells
top10.variable.genes.sr <- head(VariableFeatures(short.reads.harmony), 10)
# Convert gene IDs in top10 to their respective symbols using FeatureIDtoName() function
converted_top10.sr <- sapply(top10.variable.genes.sr, function(feature_id) {
  FeatureIDtoName(feature_id, ConvType = 'g.to.gName', gtfFilePath = gtfFilePath)
})

## Plot the 2000 most variable features and label the top10 most variable transcripts
variable.feature.plot.sr <- VariableFeaturePlot(short.reads.harmony) +
  labs(title = "Top 2000 Variable Genes across Cells\n(Integrated SR)")
variable.feature.plot.sr <- LabelPoints(plot = variable.feature.plot.sr, points = top10.variable.genes.sr, labels = converted_top10.sr,
                                        repel = TRUE, xnudge = 0, ynudge = 0)
variable.feature.plot.sr

## Create feature plots showing the expression of these genes
FoI_figures <- list()
converted_top10.sr
for (feature in top10.variable.genes.sr) {
  feature_name <- converted_top10.sr[feature]
  FoI_figures[[feature]] <- FeaturePlot(short.reads.harmony, reduction = 'umap',
                                        features = feature,
                                        order = TRUE) + ggmin::theme_powerpoint() + labs(title = feature_name)
}

grid.layout.FoI <- grid.arrange(grobs = FoI_figures, nrow = 2, ncol = 5,
                                top=textGrob("Expression Levels of Top 10 Most Variable Genes\n(Integrated Short Reads)\n",
                                             gp = gpar(fontsize = 12, fontface = "bold")))
ggsave("Integrated-SR-top10-genes.pdf", grid.layout.FoI, width = 22, height = 10)

### TOP VARIABLE ISOFORM INFORMATION - INTEGRATED LONG READS (ISO LEVEL OBJ)-------------------------
top10.isoformLvl <- head(VariableFeatures(long.reads.harmony.IsoLvl), 10)
# Convert gene IDs in top10 to their respective symbols using FeatureIDtoName() function
converted_top10.isoformlvl <- sapply(top10.isoformLvl, function(feature_id) {
  FeatureIDtoName(feature_id, ConvType = 't.to.tName', gtfFilePath = gtfFilePath, novelIsoTitle = TRUE)
})

variable.feature.plot.isoLvl <- VariableFeaturePlot(long.reads.harmony.IsoLvl) +
  labs(title = "Top 2000 Variable Isoforms across Cells")
variable.feature.plot.isoLvl <- LabelPoints(plot = variable.feature.plot.isoLvl, 
                                            points = top10.isoformLvl, labels = converted_top10.isoformlvl,
                                            repel = TRUE, xnudge = 0, ynudge = 0,
                                            max.overlaps = 100)
variable.feature.plot.isoLvl

## Compile UMAP plots showing the expression of the top10 most variable transcripts
FoI_figures <- list()
converted_top10.isoformlvl
for (feature in top10.isoformLvl) {
  feature_name <- converted_top10.isoformlvl[feature]
  FoI_figures[[feature]] <- FeaturePlot(long.reads.harmony.IsoLvl, reduction = 'umap',
                                        features = feature,
                                        order = TRUE) + ggmin::theme_powerpoint() + labs(title = feature_name)
}
grid.layout.FoI <- grid.arrange(grobs = FoI_figures, nrow = 2, ncol = 5,
                                top=textGrob("Expression Levels of Top 10 Most Variable Isoforms\n(Integrated Long Reads)\n",
                                             gp = gpar(fontsize = 12, fontface = "bold")))
ggsave("Integrated-LR-top10-isoforms.pdf", grid.layout.FoI, width = 22, height = 10)