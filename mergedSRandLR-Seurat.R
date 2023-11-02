# Script to create a merged short read and long read seurat object

setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

set.seed(4242)
library(reticulate)
use_python("/Users/manveerchuahan/miniconda3/bin/python3.11")
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)
library(grid)
library(gridExtra)





## Create a merged seurat object with both short reads and long reads then perform trajectory analysis on that------------------------
lr.geneMatrixFilePath.day25 = "/Volumes/Expansion/temp/lr-cortdiff-day25_gene_count.csv"
lr.geneMatrixFilePath.day55 = "/Volumes/Expansion/temp/lr-cortdiff-day55_gene_count.csv"

sr.matrixFilePath.day25 = "/Volumes/Expansion/CELLRANGER_counts/sr_day25_cortdiff/outs/filtered_feature_bc_matrix.h5"
sr.matrixFilePath.day55 = "/Volumes/Expansion/CELLRANGER_counts/sr_day55_cortdiff/outs/filtered_feature_bc_matrix.h5"

## Function to make sure long read count matrices are in the same format as the short read ones
prepareLRcountMatrix <- function(countmatrix.filepath){
  countmatrix <- read.csv(countmatrix.filepath, row.names = 1) %>% 
    select(-1)
  
  countmatrix$geneNames <- ListofGeneIDstoNames(rownames(countmatrix), gtfFilePath)
  
  # Move the 'geneNames' column to be the first column
  countmatrix <- countmatrix %>%
    select(geneNames, everything())
  
  # Loop through each element of the 'geneNames' column
  for (i in seq_along(countmatrix$geneNames)) {
    # If the element is NA or NULL, replace it with the corresponding row name
    if (is.null(countmatrix$geneNames[[i]])) {
      countmatrix$geneNames[[i]] <- rownames(countmatrix)[[i]]
    }
  }
  
  # Add the row sums as a new column
  countmatrix$RowSum <- rowSums(countmatrix[, -which(names(countmatrix) %in% c("geneNames"))])
  # Move the 'geneNames' column to be the first column
  countmatrix <- countmatrix %>%
    select(RowSum, everything())
  
  ## Remove any duplicate geneNames based on which has more counts
  duplicate_values <- duplicated(countmatrix$geneNames) | duplicated(countmatrix$geneNames, fromLast = TRUE)
  View(countmatrix[duplicate_values, ])
  
  # Filter to keep rows with maximum RowSum in each geneName group
  filtered_df <- countmatrix %>% 
    group_by(geneNames) %>% 
    filter(RowSum == max(RowSum)) %>% 
    slice_head(n = 1) %>% 
    ungroup()
  
  # Set row names to the geneNames column
  rownames(filtered_df) <- filtered_df$geneNames
  preserve.rownames <- rownames(filtered_df)
  # Remove the geneNames column and RowSum columns
  filtered_df$RowSum <- NULL
  filtered_df$geneNames <- NULL
  rownames(filtered_df) <- preserve.rownames
  
  filtered_df
}
formated.lr.geneMatrix.day25 <- prepareLRcountMatrix(lr.geneMatrixFilePath.day25)
formated.lr.geneMatrix.day55 <- prepareLRcountMatrix(lr.geneMatrixFilePath.day55)

## Create seurat objects for every count matrix, and make sure cells match the filtered seurat objects
## that are already being used
lr.formated.D25.seurat.obj <- CreateSeuratObject(counts = formated.lr.geneMatrix.day25, project = "LR_DAY25",
                                                 min.cells = 3, min.features = 1)
lr.formated.D55.seurat.obj <- CreateSeuratObject(counts = formated.lr.geneMatrix.day55, project = "LR_DAY55",
                                                 min.cells = 3, min.features = 1)

#*** make sure cell barcodes match the filtered seurat objects we've already made
filtered.lr.seurat.obj.D55 <- readRDS(file = "/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline/LR_DAY55-QCed.rds")
integrated.sr.seurat.obj <- readRDS(file = integrated.sr.filepath)
class(lr.formated.D55.seurat.obj)
lr.formated.D25.seurat.obj <- subset(lr.formated.D55.seurat.obj@meta.data, cells = row.names(filtered.lr.seurat.obj.D55@meta.data))
common_cells <- intersect(row.names(lr.formated.D55.seurat.obj@meta.data), row.names(filtered.lr.seurat.obj.D55@meta.data))
length(common_cells)
##