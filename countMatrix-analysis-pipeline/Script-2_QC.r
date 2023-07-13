library(Seurat)
library(tidyverse)
library(reticulate)
library(DoubletFinder)
library(gridExtra)
library(grid)

set.seed(4242)

# Initialize the variables we will be using
output_prefix = "q20_Day25"
matrixFilePath = "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/filt_q20_Day25_transcript_count.csv"
feature_type = "Isoform"

#feature_type = "Gene"
#matrixFilePath = "/data/gpfs/projects/punim0646/manveer/CELLRANGER_counts/sr_day25_cortdiff/outs/filtered_feature_bc_matrix.h5"
#matrixFilePath = "/data/gpfs/projects/punim0646/manveer/CELLRANGER_counts/sr_day55_cortdiff/outs/filtered_feature_bc_matrix.h5"
min.features = 200
max.features = 999999999
min.counts = 400
max.counts = 100000
MT_threshold = 10
npc = 10
cluster_resolution = 0.9

table1 = data.frame()
cluster_resolution.figs = list()

# Import our feature_ID_converter python script
source_python("/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/feature_ID_converter.py")
gtfFilePath <- "/data/gpfs/projects/punim0646/manveer/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

### READ IN OUR RAW COUNT FILES-----
if (tolower(feature_type) == "isoform" | tolower(feature_type) == "isoforms") {
  ## Read in the FLAMES transcript count csv file
  filt_trans_counts = matrixFilePath
  count_matrix <- read.csv(filt_trans_counts, row.names = 1)
} else if (tolower(feature_type) == "gene" | tolower(feature_type) == "genes") {
  ## Read in the cellRanger .h5 count matrix file
  count_matrix <- Read10X_h5(filename = matrixFilePath, use.names = TRUE,
                                    unique.features = TRUE)
}

### Initialize a Seurat object with raw (non-normalized) data ----

## First, we will parse features that appear in at least 3 cells, 
# and cells that contain >= 1 feature into our seurat object
seurat.obj <- CreateSeuratObject(counts = count_matrix, project = output_prefix,
                                 min.cells = 3, min.features = 1)
seurat.obj
#13600 features across 811 samples within 1 assay 

# Update table1 with number of the initial total number of cells and features
initial_cells = dim(seurat.obj)[2]
table1 <- rbind(table1, data.frame("Cells"=dim(seurat.obj)[2],
                            "Median features per cell"=median(seurat.obj$nFeature_RNA), 
                            row.names = paste0('min features > 0'),check.names = FALSE))

### Visualize metadata features to assess the quality of cells ----
# Add a column with the percentage of mitochondrial content per cell
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
#View(seurat.obj@meta.data)
median.mito.content.before = median(seurat.obj@meta.data[["percent.mt"]])

VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

## Next override our seurat object, removing unwanted cells. You can modify these parameters in the function
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > min.features & nFeature_RNA < max.features 
                     & percent.mt < MT_threshold & nCount_RNA < max.counts & nCount_RNA > min.counts)
FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

# Add the new number of cells and median features per cell to table1
table1 <- rbind(table1, data.frame("Cells"=dim(seurat.obj)[2],
                                   "Median features per cell"=median(seurat.obj$nFeature_RNA), 
                                   row.names = paste0('min features > 0', min.features),check.names = FALSE))

association.plt <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + labs(title = "Association between unique features and number of reads per cell")
#VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#vln1 <- VlnPlot(seurat.obj, features = c("nFeature_RNA")) + labs(title = "Unique Features per Cell\nAfter Removing Unwanted Cells")
#vln2 <- VlnPlot(seurat.obj, features = c("nCount_RNA")) + labs(title = "RNA Molecules per Cell\nAfter Filtering")
#vln3 <- VlnPlot(seurat.obj, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Removing Unwanted Cells")

### Now you have removed unwanted cells, it is time to normalize the data. ----
### By default : A global-scaling normalization method "LogNormalize" is  used
### that normalizes the feature expression measurements for each cell by the total expression, 
### this is multiplied by a scale factor (default = 10,000), then is log-transformed.
seurat.obj <- NormalizeData(seurat.obj)

### Next we'll find the top 2000 variable features (transcripts/genes) in our dataset ----
### Focusing our further downstream analyses on these features has been shown to be more effective
seurat.obj <- FindVariableFeatures(seurat.obj, 
                                   selection.method = 'vst',
                                   nfeatures = 2000)
# Extract the top 10 most variable features across cells
top10 <- head(VariableFeatures(seurat.obj), 10)

# Convert transcript IDs in top10 to their respective symbols using FeatureIDtoSymbol() function
#converted_top10 <- sapply(top10, function(feature_id) {
#  FeatureIDtoSymbol(feature_id, ConvType = 't.to.tName', gtfFilePath = gtfFilePath, novelIsoTitle = TRUE)
#})

if (tolower(feature_type) == "isoform" | tolower(feature_type) == "isoforms") {
  # Convert transcript IDs in top10 to their respective symbols using FeatureIDtoSymbol() function
  converted_top10 <- sapply(top10, function(feature_id) {
    FeatureIDtoSymbol(feature_id, ConvType = 't.to.tName', gtfFilePath = gtfFilePath, novelIsoTitle = TRUE)
  })
} else if (tolower(feature_type) == "gene" | tolower(feature_type) == "genes") {
  # Convert gene IDs in top10 to their respective symbols using FeatureIDtoSymbol() function
  converted_top10 <- sapply(top10, function(feature_id) {
    FeatureIDtoSymbol(feature_id, ConvType = 'g.to.gName', gtfFilePath = gtfFilePath)
  })
} else {
  converted_top10 <- top10
}

# Plot the 2000 most variable features and label the top10 most variable transcripts
variable.feature.plot <- VariableFeaturePlot(seurat.obj) +
  labs(title = "Top 2000 Variable Features across Cells")
variable.feature.plot <- LabelPoints(plot = variable.feature.plot, points = top10, labels = converted_top10,
            repel = TRUE, xnudge = 0, ynudge = 0)
variable.feature.plot

### Next we apply a linear transformation (scaling) ----
### this is a standard pre-processing step prior to dimensional reduction
all.features <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.features)
### Next apply linear dimensional reduction with PCA, using our variable features ----
seurat.obj <- RunPCA(seurat.obj,
                     features = VariableFeatures(object = seurat.obj))
# Use an elbow plot to determine the dimensionality of our data
elbow.plt <- ElbowPlot(seurat.obj) + labs(title = 'Standard Deviation explained by each Principle Component')
### Now we'll cluster our cells ----

# FindNeightbors() makes a KNN graph based on euclidean distance in PCA space and 
# refines edge weights between any two cells based on the shared overlap in their 
# local neighborhoods (jaccard similarity). 
# We'll be using dimensions based on the elbowplot just generated as the input
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npc)

# Now we actually cluster our cells, FindClusters() needs a resolution parameter 
# the resolution sets the granularity of downstream clustering. 
# Settings are recommended between 0.4-1.2, but may need to be increased for larger datasets

seurat.obj <- FindClusters(seurat.obj, resolution = c(0.3, 0.5, 0.7, 0.9, 1.1, cluster_resolution))

# save each plot to review in the final pdf
cluster_resolution.figs[[1]] <- DimPlot(seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE)
cluster_resolution.figs[[2]] <- DimPlot(seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
cluster_resolution.figs[[3]] <- DimPlot(seurat.obj, group.by = "RNA_snn_res.0.7", label = TRUE)
cluster_resolution.figs[[4]] <- DimPlot(seurat.obj, group.by = "RNA_snn_res.0.9", label = TRUE)
cluster_resolution.figs[[5]] <- DimPlot(seurat.obj, group.by = "RNA_snn_res.1.1", label = TRUE)
desired_res_colTitle <- paste0("RNA_snn_res.", cluster_resolution)
cluster_resolution.figs[[6]] <- DimPlot(seurat.obj, group.by = desired_res_colTitle, label = TRUE) + 
  labs(title = paste0("Desired cluster resolution : ", cluster_resolution))


# Make sure we've set the correct cluster identities based on our desired cluster_resolution
desired_resolution <- paste0("RNA_snn_res.", cluster_resolution)
Idents(seurat.obj) <- desired_resolution
#Idents(seurat.obj)

### Move on the non-linear dimensionality reduction (UMAP)----
seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc)


### filter out doublets (remember  to modify doublet rate if samples have variable target cells)

### Next move on to Doublet detection with DoubletFinder -----
## First identify the most optimum pK value (no-ground truth strategy)
sweep.res.list <- paramSweep_v3(seurat.obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
#View(sweep.stats)
BCmvn <- find.pK(sweep.stats)

pK <- BCmvn %>% # select the pK that corresponds to max BCmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 
annotations <- seurat.obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*nrow(seurat.obj@meta.data))  ## Assuming 3.9% doublet formation rate - based on original pool of 5000 cells
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Then run DoubletFinder 
seurat.obj <- doubletFinder_v3(seurat.obj, 
                                  PCs = 1:20, 
                                  pN = 0.25, 
                                  pK = pK, 
                                  nExp = nExp_poi.adj,
                                  reuse.pANN = FALSE, sct = FALSE)

# Clean up DoubletFinder's classification column in seurat object's metadata
colnames(seurat.obj@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(seurat.obj@meta.data))
# Create a summary table showing doublet information from our dataset
statsDoublets <- seurat.obj@meta.data %>% 
  group_by(DF.classifications) %>%
  summarize(median_nCount_RNA = median(nCount_RNA),
            median_nFeature_RNA = median(nFeature_RNA),
            count = n())

# Visualize doublets in our UMAP plot
doublets.umap <- DimPlot(seurat.obj, reduction = 'umap', group.by = "DF.classifications") + 
  labs(title = "Doublets that were detected and removed") +
  ggmin::theme_powerpoint()

# Only keep the cells that are considered singlets - discard the doublets
seurat.obj <- subset(seurat.obj, subset = DF.classifications == 'Singlet')

### Generate QC Violin Plots for each cluster----
vln1 <- VlnPlot(seurat.obj, features = c("nFeature_RNA")) + labs(title = "Unique Features per Cell\nAfter Removing Unwanted Cells") +
  geom_hline(yintercept = median(seurat.obj$nFeature_RNA), linetype = 'dashed')
vln2 <- VlnPlot(seurat.obj, features = c("nCount_RNA")) + labs(title = "RNA Molecules per Cell\nAfter Filtering") +
  geom_hline(yintercept = median(seurat.obj$nCount_RNA), linetype = 'dashed')
vln3 <- VlnPlot(seurat.obj, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Removing Unwanted Cells") +
  geom_hline(yintercept = median(seurat.obj@meta.data[["percent.mt"]]), linetype = 'dashed')
#association.plt <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
#  geom_smooth(method = "lm")

### Plot the final UMAP and save some statistics on the final filtered seurat object----
umap_title <- paste0(output_prefix, " Sample\n(Resolution : ", cluster_resolution, ", Dimensions: ", npc)
baseUMAPplot <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, ) + 
  labs(color = "cluster \n(from PCA)", 
       title = umap_title) +
  ggmin::theme_powerpoint()

filtered.summary <- rbind("Sample ID" = output_prefix,
                          "Initial Cell Count" = initial_cells,
                          "Cells after QC" = dim(seurat.obj)[2],
                          "Feature Type" = feature_type,
                          "Min. Feature Threshold" = min.features,
                          "Max. Feature Threshold" = max.features,
                          "Min. Counts Threshold" = min.counts,
                          "Max. Counts Threshold" = max.counts,
                          "Median Features per Cell" = median(seurat.obj$nFeature_RNA),
                          "Median Reads per Cell"= median(seurat.obj$nCount_RNA),
                          "Mitochondrial Content Threshold (%)" = MT_threshold,
                          "Median MT % before filter" = median.mito.content.before,
                          "Median MT % after filter" = median(seurat.obj@meta.data[["percent.mt"]]),
                          "NPCs (UMAP dimensions)" = npc)
                         

### Compile and arrange all the figures that were generated ----
baseUMAPplot
elbow.plt
filtered.summary
association.plt
vln1 | vln2 | vln3
doublets.umap
statsDoublets

final.layout.alt <- grid.arrange(tableGrob(filtered.summary), baseUMAPplot,
                               elbow.plt, variable.feature.plot,
                               vln1, vln2,
                               vln3, association.plt,
                               doublets.umap, tableGrob(statsDoublets),
                               nrow = 5, ncol = 2,
                               top = textGrob(paste0("SCRIPT 2 (QC) SUMMARY (", output_prefix,")\n"),
                                              gp = gpar(fontsize = 20)))

pdf.figs.title <- paste0(output_prefix, "-QC-summary-script2.pdf")
ggsave(pdf.figs.title, final.layout.alt, width = 20, height = 30)

cluster.res.layout <- grid.arrange(cluster_resolution.figs[[1]], cluster_resolution.figs[[2]],
                                   cluster_resolution.figs[[3]], cluster_resolution.figs[[4]], 
                                   cluster_resolution.figs[[5]], cluster_resolution.figs[[6]], 
                                   tableGrob(filtered.summary),
                                   nrow = 3, ncol = 3,
                                   top=textGrob(paste0("SCRIPT 2 (QC) - ", output_prefix, " CLUSTER RESOLUTION FIGURES"),
                                                gp = gpar(fontsize = 15)))
res.figs.title <- paste0(output_prefix, "-cluster-res-figures-script2.pdf")
ggsave(res.figs.title, cluster.res.layout, width = 20, height = 18)


### Compile UMAP plots showing the expression of the top10 most variable features----
FoI_figures <- list()
converted_top10
for (feature in top10) {
  feature_name <- converted_top10[feature]
  FoI_figures[[feature]] <- FeaturePlot(seurat.obj, reduction = 'umap',
                                        features = feature,
                                        order = TRUE) + ggmin::theme_powerpoint() + labs(title = feature_name)
}
grid.layout.FoI <- grid.arrange(grobs = FoI_figures, nrow = 2, ncol = 5,
                                top=textGrob(paste0("Expression Levels of Top 10 Most Variable Features (", output_prefix,")\n"),
                                gp = gpar(fontsize = 12, fontface = "bold")))
figs.title <- paste0(output_prefix, "-top10-features-script2.pdf")
ggsave(figs.title, grid.layout.FoI, width = 22, height = 10)


### Export the processed seurat object for use in further scripts
rdsFilename <- paste0(output_prefix,"-QCed.rds")
rdsFilename
saveRDS(seurat.obj, file = rdsFilename)
