library(celldex)
library(SingleR)
library(Seurat)
library(tidyverse)
library(pheatmap)


# How to use sctype in R : https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md


# Get reference dataset from the celldex package-----
reference <- celldex::HumanPrimaryCellAtlasData()
#View(as.data.frame(colData(reference)))

# Load the Seurat object (from script 2) using its RDS file-----
seurat.obj <- readRDS(file = "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/filtered_sr_Day25-QCed.rds")

### run SingleR (default mode)-------------------------------------------
# default for SingleR is to perform annotation of each individual cell in the test data

# The counts slot contains the raw counts used in your seurat object
raw_counts <- GetAssayData(seurat.obj, slot = 'counts')

prediction <- SingleR(test = raw_counts,
                      ref = reference,
                      labels = reference$label.fine)
prediction

# Add the singleR predictions to the seurat object metadata, making sure labels are matched to their respective barcodes
seurat.obj$singleR.labels <- prediction$labels[match(rownames(seurat.obj@meta.data), rownames(prediction))]
DimPlot(seurat.obj, reduction = "umap", group.by = "singleR.labels")

### Annotation Diagnostics -------------------------------------------

## Based on the scores in singleR output for each cell
prediction$scores
plotScoreHeatmap(prediction)

## Based on deltas across cells
plotDeltaDistribution(prediction)

## Comparing to unsupervised clustering

# create a table calculating the numbers of cell type labels for each cluster
tab <- table(Assigned = prediction$labels, Clusters = pbmc.seurat.filtered$seurat_clusters)

# The lower values should be white, the higher values should be blue and there should be
# 10 colour values between those extremes
pheatmap(log10(tab + 10), color = colorRampPalette(c('white', 'blue'))(10))

