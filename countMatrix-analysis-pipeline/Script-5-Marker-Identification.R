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

## Read in filtered seurat objects for both samples with short reads and long reads---------------
integrated.lr.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/Relabeled_LR_days25+55-script3.rds"
integrated.sr.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/Relabeled_SR_days25+55-script3.rds"
integrated.lr.seurat.obj <- readRDS(file = integrated.lr.filepath)
integrated.sr.seurat.obj <- readRDS(file = integrated.sr.filepath)

## REMOVE COLS IN METADATA THAT AREN'T NEEDED ANYMORE------------------------
integrated.lr.seurat.obj@meta.data$RNA_snn_res.0.3 <- NULL
integrated.lr.seurat.obj@meta.data$RNA_snn_res.0.5 <- NULL
integrated.lr.seurat.obj@meta.data$RNA_snn_res.0.7 <- NULL
integrated.lr.seurat.obj@meta.data$RNA_snn_res.0.9 <- NULL
integrated.lr.seurat.obj@meta.data$RNA_snn_res.1.1 <- NULL

# Create function to get conserved markers for any given cluster----------
get_conserved.LR <- function(cluster){
  FindConservedMarkers(integrated.lr.seurat.obj,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
get_cluster_top10_GeneMarkers.LR <- function(cluster){
  conserved.markers <- map_dfr(cluster, get_conserved.LR)
  
  top10 <- conserved.markers %>% 
    mutate(avg_fc = (LR_DAY25_avg_log2FC + LR_DAY55_avg_log2FC) /2) %>% 
    group_by(cluster_id) %>% 
    top_n(n = 15,
          wt = avg_fc)
  
  correspondingGeneNames <- ListofGeneIDstoNames(top10$gene, gtfFilePath)
  top10$geneName <- correspondingGeneNames
  
  top10
}
generate.DimPlots.FOI.LR <- function(seurat.obj, FOIs, cluster){
  FoI_figures <- list()
  rows <- ceiling(length(FOIs)/2)
  
  # Convert gene IDs in Features of Interest (FOIs) list to their respective symbols using ListofGeneIDstoNames() function
  converted_geneNames <- ListofGeneIDstoNames(FOIs, gtfFilePath = gtfFilePath)

  for (index in seq_along(FOIs)) {
    feature_name <- converted_geneNames[index]
    FoI_figures[[index]] <- FeaturePlot(object = seurat.obj, 
                features = FOIs[[index]],
                order = TRUE,
                min.cutoff = 'q10', 
                label = TRUE,
                repel = TRUE) +
      labs(title = feature_name)
  }
  
  grid.layout.FoI <- grid.arrange(grobs = FoI_figures, nrow = rows, ncol = 2,
                                  top=textGrob(paste0("Cluster Markers of Interest - Cluster ", cluster,"\n"),
                                               gp = gpar(fontsize = 12, fontface = "bold")))
  
  grid.layout.FoI
}

## Identification of conserved markers in all conditions LONG READ data-------------------------
DimPlot(integrated.lr.seurat.obj, reduction = 'umap', group.by = 'customclassif',
        label = TRUE)
DimPlot(integrated.lr.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters',
        label = TRUE)
DimPlot(integrated.lr.seurat.obj, reduction = 'umap', group.by = 'orig.ident')

## INVESTIGATE CLUSTER 8 - DGE--------------------------------------
# Extract top 10 markers per cluster
cluster8_top10 <- get_cluster_top10_GeneMarkers.LR(8)

View(cluster8_top10)

# Plot interesting marker gene expression for cluster 8 (SLC17A6, STMN2, STMN4, GAP43, INA)
cluster8_FoIs <- list("ENSG00000091664.9", "ENSG00000104435.14", 
                   "ENSG00000015592.17", "ENSG00000172020.13", 
                   "ENSG00000148798.11")

generate.DimPlots.FOI.LR(integrated.lr.seurat.obj, cluster8_FoIs, 8)

FeaturePlot(object = integrated.lr.seurat.obj, 
            features = c("ENSG00000091664.9", "ENSG00000104435.14", 
                         "ENSG00000015592.17", "ENSG00000172020.13", 
                         "ENSG00000148798.11"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE,
            split.by = 'orig.ident')

VlnPlot(object = integrated.lr.seurat.obj, 
        features = c("ENSG00000091664.9", "ENSG00000104435.14", 
                     "ENSG00000015592.17", "ENSG00000172020.13", 
                     "ENSG00000148798.11"))

## INVESTIGATE CLUSTERS 9 AND 6 (hypothesized to be neural stem cells)-----------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$seurat_clusters
neuralstemcells.top10 <- get_cluster_top10_GeneMarkers.LR(c(9,6))
View(neuralstemcells.top10)

## WHATS DIFFERENTIALLY EXPRESSED BETWEEN CLUSTER 8 AND CLUSTER 5