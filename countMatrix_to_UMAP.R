# SCRIPT THAT GENERATES UMAP (and other) PLOTS FOR BLAZE+FLAMES GENERATED COUNT MATRICES

# Load libraries
#BiocManager::install("Seurat")
#reticulate::py_install(packages = 'umap-learn')
#install.packages("tidyverse")
library(Seurat)
library(tidyverse)
install.packages("gridExtra")

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
# thus 0.7 seems appropriate for this Q20 dataset

# setting the resolution for cell cluster identities
Idents(seurat.obj) <- "RNA_snn_res.0.7"
Idents(seurat.obj)

## STEP 8 : NON-LINEAR DIMENSIONALITY REDUCTION --------------------------------
seurat.obj <- RunUMAP(seurat.obj, dims = 1:10)
baseUMAPplot <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, ) + 
  labs(color = "cluster \n(from PCA)", 
       title = 'Q20 Day 25 cort-diff Sample \n(Resolution: 0.7, Dimensions: 10')

#Install an alternative (cleaner) Seurat-applied theme called ggmin
#install.packages("devtools")  
#devtools::install_github("sjessa/ggmin")
baseplot + ggmin::theme_powerpoint()

## STEP 9 : FIND CLUSTER MARKERS------------------------------------------------
q20_cluster_markers <- FindAllMarkers(seurat.obj, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)
q20_cluster_markers

## STEP 10 : GENERATE TRANSCRIPT EXPRESSION UMAPS-------------------------------
UMIcountUMAP <- FeaturePlot(seurat.obj, reduction = "umap", features = 'nCount_RNA')+
  labs(color = "UMI count",title = '')+ theme(text = element_text(size = 10))
UMIcountUMAP + ggmin::theme_powerpoint()

FeatureCountUMAP <- FeaturePlot(seurat.obj, reduction = "umap", features = 'nFeature_RNA')+
  labs(color = str_wrap("Isoform count",15),title = '')+ theme(text = element_text(size = 10))
FeatureCountUMAP + ggmin::theme_powerpoint()

top10
isoforms_of_interest = c("ENSG00000104435.14-79611117-79665011-1", 
                         "ENST00000314744.6", 
                         "ENSG00000131711.15-72107336-72195422-1")
isoforms_of_interest_figures <- list()
for (isoform in isoforms_of_interest) {
  isoforms_of_interest_figures <- append(rst_figures, 
                                         FeaturePlot(seurat.obj, reduction = 'umap',
                                                                 features = isoform,
                                                                 order = TRUE) + ggmin::theme_powerpoint())
}
grid.
