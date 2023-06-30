# SCRIPT THAT GENERATES UMAP (and other) PLOTS FOR BLAZE+FLAMES GENERATED COUNT MATRICES

# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)

# Create a list to store the plots that we generate
global_figures = list()
additional_figures = list()

### STEP 0 : CREATING SEURAT OBJECT----------------------------------------------
## Read in the transcript count csv file
count_matrix <- read.csv("/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_1000expCells/transcript_count.csv",
                         row.names = 1)
# Remove gene_id column from count matrix
count_matrix <- subset(count_matrix, select = -gene_id)

transcript_gene_key <- read.csv("/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_1000expCells/transcript_count.csv")
transcript_gene_key <- transcript_gene_key[,c("transcript_id", "gene_id")]

## Initialize a Seurat object with raw (non-normalized) data
# We are only parsing features that appear in at least 3 cells, 
# and cells that contain >= 200 features into our seurat object
seurat.obj <- CreateSeuratObject(counts = count_matrix, project = "Q20",
                                 min.cells = 3, min.features = 200)
str(seurat.obj)
seurat.obj
# 16772 features across 796 samples within 1 assay

### STEP 1 : QUALITY CONTROL-----------------------------------------------------
# View seurat object's metadata slot
View(seurat.obj@meta.data)

# Visualize metadata features to assess the quality of cells
VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

### STEP 2 : FILTERING-----------------------------------------------------------
# all cells seemed good quality based on the previous QC step
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 999999) 
seurat.obj

### STEP 3 : NORMALIZE DATA------------------------------------------------------
seurat.obj <- NormalizeData(seurat.obj)
str(seurat.obj)

### STEP 4 : IDENTIFY HIGHLY VARIABLE FEATURES ----------------------------------
seurat.obj <- FindVariableFeatures(seurat.obj, 
                                   selection.method = 'vst',
                                   nfeatures = 2000)
# Identify the top 10 most highly variable transcripts across cells
top10 <- head(VariableFeatures(seurat.obj), 10)
top10

# Plot the 2000 most variable features and label the top10 most variable transcripts
plot1 <- VariableFeaturePlot(seurat.obj)
LabelPoints(plot = plot1, points = top10)


### STEP 5 : SCALING-------------------------------------------------------------
all.features <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.features)

### STEP 6 : LINEAR DIMENSIONALITY REDUCTION-------------------------------------
seurat.obj <- RunPCA(seurat.obj,
                     features = VariableFeatures(object = seurat.obj))
# Look at the first 5 Principle components and the top 5 features in each PC
print(seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

# Determine the dimensionality of our data
ElbowPlot(seurat.obj) # It seems that 15 PCs are sufficient to explain the variation within our data
global_figures[[1]] <- ElbowPlot(seurat.obj)

### STEP 7 : CLUSTERING ---------------------------------------------------------
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:10)

# Determine the correct resolution for our clusters
seurat.obj <- FindClusters(seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
View(seurat.obj@meta.data)

DimPlot(seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(seurat.obj, group.by = "RNA_snn_res.0.7", label = TRUE)
DimPlot(seurat.obj, group.by = "RNA_snn_res.1", label = TRUE)

# With a resolution of 0.7, cells in each clusters seem distinctly separated
# thus 0.7 seems appropriate for this Q20 dataset. ****ADD 0.7 PLOT to the list*****

# setting the resolution for cell cluster identities
Idents(seurat.obj) <- "RNA_snn_res.0.7"
Idents(seurat.obj)

### STEP 8 : NON-LINEAR DIMENSIONALITY REDUCTION --------------------------------
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
baseUMAPplot <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, ) + 
  labs(color = "cluster \n(from PCA)", 
       title = 'Q20 Day 25 cort-diff Sample \n(Resolution: 0.7, Dimensions: 10') +
  ggmin::theme_powerpoint()
#baseUMAPplot
global_figures[[2]] <- baseUMAPplot

### STEP 9 : FIND CLUSTER MARKERS------------------------------------------------
q20_cluster_markers <- FindAllMarkers(seurat.obj, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)
q20_cluster_markers

### STEP 10 : GENERATE TRANSCRIPT EXPRESSION UMAPS-------------------------------

## Create a UMAP plot showing the number of RNA molecules per cell
UMIcountUMAP <- FeaturePlot(seurat.obj, reduction = "umap", features = 'nCount_RNA')+
  labs(color = "UMI count",title = '')+ theme(text = element_text(size = 10)) + ggmin::theme_powerpoint()
#UMIcountUMAP
global_figures[[3]] <- UMIcountUMAP
#global_figures[[3]]


## Create a UMAP plot showing the number of features expressed per cell
FeatureCountUMAP <- FeaturePlot(seurat.obj, reduction = "umap", features = 'nFeature_RNA')+
  labs(color = str_wrap("Isoform count",15),title = '')+ theme(text = element_text(size = 10)) +
  ggmin::theme_powerpoint()
#FeatureCountUMAP 
global_figures[[4]] <- FeatureCountUMAP

## Create a UMAP plot showing the expression levels for specified isoforms of interest
top10
# Here I chose three different isoforms from the top10 most variable isoform list (commented above)
isoforms_of_interest = c("ENSG00000104435.14-79611117-79665011-1", 
                         "ENST00000314744.6", 
                         "ENSG00000131711.15-72107336-72195422-1",
                         "ENSG00000149591.17-117199370-117204778-1")

# Genes that correspond to isoforms listed above : 
# STMN2 [2], NEUROG1 [8], MAP1B [9], TAGLN [1]
IoI_figures <- list()
for (isoform in isoforms_of_interest) {
  IoI_figures[[isoform]] <- FeaturePlot(seurat.obj, reduction = 'umap', 
                                        features = isoform,
                                        order = TRUE) + ggmin::theme_min()
}

# Create the grid layout to view expression levels of the Isoforms of interest
grid_layout_IoI <- grid.arrange(grobs = IoI_figures, ncol = 2, 
                                top=textGrob("Isoform Expression levels across Cells \n(Taken from Top 10 most variable features)\n", 
                                             gp = gpar(fontsize = 12, fontface = "bold")))

             
### STEP 11 : COMPILE PLOTS IN GRID FORMAT--------------------------------------

## Create a grid layout for the global_figures list
grid_layout_main <- grid.arrange(grobs = global_figures, ncol = 2)

## Concatenate layout with isoforms of interest plots
concatenated_layout <- grid.arrange(grid_layout_IoI, grid_layout_main)
