setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

# import libraries
library(tidyverse)
library(reticulate)
library(png)
library(gridExtra)
library(grid)
library(Seurat)
library(UpSetR)

# Set seed for reproducibility
set.seed(4242)

# Import our python script that creates heatmap with seaborn
use_python("/Users/manveerchuahan/miniconda3/bin/python3.11")
source_python("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline/seaborn_matrix.py")

## Read in filtered + integrated seurat objects for both samples with short reads and long reads---------------
integrated.lr.filepath <- "Relabeled_LR_days25+55-script3.rds"
integrated.sr.filepath <- "Relabeled_SR_days25+55-script3.rds"
integrated.sr.regressed.filepath <- "INTEGRATED_SR_days25_55-REGRESSED.rds"
integrated.lr.regressed.filepath <- "INTEGRATED_LR_days25+55-REGRESSED.rds"


integrated.lr.seurat.obj <- readRDS(file = integrated.lr.filepath)
integrated.sr.seurat.obj <- readRDS(file = integrated.sr.filepath)
integrated.lr.seurat.obj.regressed <- readRDS(file = integrated.lr.regressed.filepath)
integrated.sr.seurat.obj.regressed <- readRDS(file = integrated.sr.regressed.filepath)

## UMAPS WITH LABELS-------------------------
lr.celltype <- DimPlot(integrated.lr.seurat.obj, reduction = 'umap', group.by = 'customclassif',
                       label = TRUE, label.size = 4, repel = TRUE) +
  labs(title = 'Integrated Long Reads\nscType Generated Labels')
sr.celltype <- DimPlot(integrated.sr.seurat.obj, reduction = 'umap', group.by = 'customclassif',
                       label = TRUE, label.size = 4, repel = TRUE) +
  labs(title = 'Integrated Short Reads\nscType Generated Labels')

lr.unsupervised <- DimPlot(integrated.lr.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters',
                           label = TRUE, label.size = 5, repel = TRUE) +
  labs(title = 'Integrated Long Reads\nUnsupervised Seurat Clustering')
sr.unsupervised <- DimPlot(integrated.sr.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters',
                           label = TRUE, label.size = 5, repel = TRUE) +
  labs(title = 'Integrated Short Reads\nUnsupervised Seurat Clustering')

DimPlot(integrated.lr.seurat.obj, reduction = 'umap', group.by = 'orig.ident')

sr.celltype | sr.unsupervised

lr.celltype | lr.unsupervised


#lr.celltype.regressed <- DimPlot(integrated.lr.seurat.obj.regressed, reduction = 'umap', group.by = 'customclassif',
#                       label = TRUE, label.size = 4, repel = TRUE) +
#  labs(title = 'Integrated Long Reads\nscType Generated Labels')
#sr.celltype.regressed <- DimPlot(integrated.sr.seurat.obj.regressed, reduction = 'umap', group.by = 'customclassif',
#                       label = TRUE, label.size = 4, repel = TRUE) +
#  labs(title = 'Integrated Short Reads\nscType Generated Labels')

#lr.unsupervised.regressed <- DimPlot(integrated.lr.seurat.obj.regressed, reduction = 'umap', group.by = 'seurat_clusters',
#                           label = TRUE, label.size = 5, repel = TRUE) +
#  labs(title = 'Integrated Long Reads\nUnsupervised Seurat Clustering')
#sr.unsupervised.regressed <- DimPlot(integrated.sr.seurat.obj.regressed, reduction = 'umap', group.by = 'seurat_clusters',
#                           label = TRUE, label.size = 5, repel = TRUE) +
#  labs(title = 'Integrated Short Reads\nUnsupervised Seurat Clustering')

#DimPlot(integrated.lr.seurat.obj.regressed, reduction = 'umap', group.by = 'orig.ident')

#sr.unsupervised.regressed | lr.unsupervised.regressed

### Define functions to calculate shared cells and produce heatmaps-------------------------------
calculateSharedCellsPerCluster <- function(LR.seuratObj, SR.seuratObj){
  # First identify barcodes that are only found in LR and absent in SR object (to remove later)
  LR.barcodes <- rownames(LR.seuratObj@meta.data) %>% substr(1, nchar(.) - 2)
  SR.barcodes <- rownames(SR.seuratObj@meta.data) %>% substr(1, nchar(.) - 2)
  LR.exclusive.barcodes <- setdiff(LR.barcodes, SR.barcodes)
  print(length(LR.exclusive.barcodes))
  
  num.unique.Clusters.LR <- length(unique(LR.seuratObj@meta.data$seurat_clusters))
  num.unique.Clusters.SR <- length(unique(SR.seuratObj@meta.data$seurat_clusters))
  
  # Initialize an empty matrix to store the results
  overlap_matrix <- matrix(0, nrow = num.unique.Clusters.LR, ncol = num.unique.Clusters.SR)
  colnames(overlap_matrix) <- 0:(num.unique.Clusters.SR - 1)
  rownames(overlap_matrix) <- 0:(num.unique.Clusters.LR - 1)

  
  for(LRcluster in 0:(num.unique.Clusters.LR - 1)){
    cluster_subset_LR <- LR.seuratObj@meta.data %>% filter(seurat_clusters == LRcluster)
    # Remove the last two digits from row names
    rownames(cluster_subset_LR) <- substr(rownames(cluster_subset_LR), 1, nchar(rownames(cluster_subset_LR)) - 2)
    # Remove barcodes that are exclusive to long reads and not found in short reads 
    # as (this skews our calculate proportions)
    cluster_subset_LR <- cluster_subset_LR[!(rownames(cluster_subset_LR) %in% LR.exclusive.barcodes), ]
    
    for(SRcluster in 0:(num.unique.Clusters.SR - 1)){
      cluster_subset_SR <- SR.seuratObj@meta.data %>% filter(seurat_clusters == SRcluster)
      # Remove the last two digits from row names
      rownames(cluster_subset_SR) <- substr(rownames(cluster_subset_SR), 1, nchar(rownames(cluster_subset_SR)) - 2)
      
      
      # Calculate overlapping cells
      overlapping_cells <- length(intersect(rownames(cluster_subset_LR), rownames(cluster_subset_SR)))
      
      # Calculate proportion of long read cells found in short read object
      if(nrow(cluster_subset_LR) != 0) {
        prop_overlap_LR <- round(overlapping_cells / nrow(cluster_subset_LR) * 100)
      } else {
        prop_overlap_LR <- NA
      }
      
      # Store results in the matrix
      overlap_matrix[LRcluster + 1, SRcluster + 1] <- prop_overlap_LR
    }
  }

  return(overlap_matrix)
}

calculateSharedCellsPerLabel <- function(LR.seuratObj, SR.seuratObj){
  # First identify barcodes that are only found in LR and absent in SR object (to remove later)
  LR.barcodes <- rownames(LR.seuratObj@meta.data) %>% substr(1, nchar(.) - 2)
  SR.barcodes <- rownames(SR.seuratObj@meta.data) %>% substr(1, nchar(.) - 2)
  LR.exclusive.barcodes <- setdiff(LR.barcodes, SR.barcodes)
  print(length(LR.exclusive.barcodes))
  
  num.unique.Clusters.LR <- length(unique(LR.seuratObj@meta.data$customclassif))
  LR.cellLabel.list <- unique(LR.seuratObj@meta.data$customclassif)
  
  num.unique.Clusters.SR <- length(unique(SR.seuratObj@meta.data$customclassif))
  SR.cellLabel.list <- unique(SR.seuratObj@meta.data$customclassif)
  
  # Initialize an empty matrix to store the results
  overlap_matrix <- matrix(0, nrow = num.unique.Clusters.LR, ncol = num.unique.Clusters.SR)
  colnames(overlap_matrix) <- SR.cellLabel.list
  rownames(overlap_matrix) <- LR.cellLabel.list
  
  for(LRcluster in LR.cellLabel.list){
    cluster_subset_LR <- LR.seuratObj@meta.data %>% filter(customclassif == LRcluster)
    # Remove the last two digits from row names
    rownames(cluster_subset_LR) <- substr(rownames(cluster_subset_LR), 1, nchar(rownames(cluster_subset_LR)) - 2)
    # Remove barcodes that are exclusive to long reads and not found in short reads 
    # as (this skews our calculate proportions)
    cluster_subset_LR <- cluster_subset_LR[!(rownames(cluster_subset_LR) %in% LR.exclusive.barcodes), ]
    
    for(SRcluster in SR.cellLabel.list){
      cluster_subset_SR <- SR.seuratObj@meta.data %>% filter(customclassif == SRcluster)
      # Remove the last two digits from row names
      rownames(cluster_subset_SR) <- substr(rownames(cluster_subset_SR), 1, nchar(rownames(cluster_subset_SR)) - 2)
      
      
      # Calculate overlapping cells
      overlapping_cells <- length(intersect(rownames(cluster_subset_LR), rownames(cluster_subset_SR)))
      
      # Calculate proportion of long read cells found in short read object
      if(nrow(cluster_subset_LR) != 0) {
        prop_overlap_LR <- round(overlapping_cells / nrow(cluster_subset_LR) * 100)
      } else {
        prop_overlap_LR <- NA
      }
      
      # Store results in the matrix
      overlap_matrix[LRcluster, SRcluster] <- prop_overlap_LR
    }
  }
  
  return(overlap_matrix)
}

## Generate a heatmap showing the proportions of shared cells between unsupervised clusters in LR and SR
results2 <- calculateSharedCellsPerCluster(integrated.lr.seurat.obj, integrated.sr.seurat.obj)
plotHeatmapfrom2DMatrix(results2, filename = "cluster-overlap-heatmap-noRegression.png", 
                        plotTitle = "Percent overlap between LR and SR Unsupervised Clusters")
img <- readPNG("cluster-overlap-heatmap-noRegression.png")
grid.newpage()
grid::grid.raster(img)

sr.unsupervised | lr.unsupervised


## Generate heatmap displaying the shared cells between cell type labels (from scType) in LR and SR
results3 <- calculateSharedCellsPerLabel(integrated.lr.seurat.obj, integrated.sr.seurat.obj)
plotHeatmapfrom2DMatrix(results3, filename = "cell-label-overlap-heatmap-noRegression.png",
                        plotTitle = "Percent overlap between labeled Cell Types in LR and SR",
                        useCellLabels = TRUE,
                        LRcellTypeLabels = unique(integrated.lr.seurat.obj$customclassif),
                        SRcellTypeLabels = unique(integrated.sr.seurat.obj@meta.data$customclassif),
                        labelTilt = TRUE,
                        biggerPlotSize = TRUE,
                        colour = "Reds")
img2 <- readPNG("cell-label-overlap-heatmap-noRegression.png")
grid.newpage()
grid::grid.raster(img2)
sr.celltype | lr.celltype

### Make an upset plot-----------------------------------
# Create a data frame for LR
lr_df <- data.frame(
  Cell = colnames(integrated.lr.seurat.obj),
  Cluster = integrated.lr.seurat.obj@meta.data$customclassif
)
# Remove last two characters from every cell barcode
lr_df$Cell <- substr(lr_df$Cell, 1, nchar(lr_df$Cell) - 2)

# Create a data frame for SR
sr_df <- data.frame(
  Cell = colnames(integrated.sr.seurat.obj),
  Cluster = integrated.sr.seurat.obj@meta.data$customclassif
)
# Remove last two characters from every cell barcode
sr_df$Cell <- substr(sr_df$Cell, 1, nchar(sr_df$Cell) - 2)

View(lr_df)
View(sr_df)

lr_df$Cluster <- paste0("LR_", lr_df$Cluster)
sr_df$Cluster <- paste0("SR_", sr_df$Cluster)

both <- rbind(lr_df,sr_df)

both_list <- split(both, both$Cluster)
View(both_list)
both_list2 <- list(`LR Neural Progenitor Cells`=both_list$`LR_Neural Progenitor cells`$Cell,
                   `SR Neural Progenitor Cells`=both_list$`SR_Neural Progenitor cells`$Cell,
                   
                   `LR Endothelial Cells`=both_list$`LR_Endothelial cells`$Cell,
                   `SR Endothelial Cells`=both_list$`SR_Endothelial cells`$Cell,
                   
                   `LR Immature Neurons`=both_list$`LR_Immature neurons`$Cell,
                   `SR Immature Neurons`=both_list$`SR_Immature neurons`$Cell,
                   
                   `LR Radial Glial Cells`=both_list$`LR_Radial glial cells`$Cell,
                   `SR Radial Glial Cells`=both_list$`SR_Radial glial cells`$Cell,
                   
                   `LR Mature Neurons`=both_list$`LR_Mature neurons`$Cell,
                   `SR Mature Neurons`=both_list$`SR_Mature neurons`$Cell,

                   `SR Cancer Cells`=both_list$`SR_Cancer cells`$Cell,
                   `SR GABAergic Neurons`=both_list$`SR_GABAergic neurons`$Cell,
                   `LR Myelinating Schwann Cells`=both_list$`LR_Myelinating Schwann cells`$Cell,
                   `SR Neuroblasts`=both_list$SR_Neuroblasts$Cell)

upset.plt <- UpSetR::upset(fromList(both_list2), nsets = length(both_list2))
upset.plt

## Export figures as a pdf
pdf("Cluster-overlap-Script6-p1.pdf", width = 11, height = 8)
grid::grid.raster(img)
grid::grid.raster(img2)
upset.plt
dev.off()

pdf("Cluster-overlap-Script6-p2.pdf", width = 11, height = 8)
grid.draw(sr.unsupervised | lr.unsupervised)
grid.draw(sr.celltype | lr.celltype)
dev.off()