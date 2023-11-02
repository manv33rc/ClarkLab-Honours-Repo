setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)
library(reticulate)
library(ComplexHeatmap)
library(circlize)
library(stringr)

# Set seed for reproducibility
set.seed(4242)

# Import our feature_ID_converter python script
source_python("feature_ID_converter.py")
gtfFilePath <- "/Volumes/Expansion/temp/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

### GENERATE A TRANSCRIPT - GENE "DICTIONARY" DERIVED FROM THE OUTPUT OF THE BLAZE PAPER'S FLAMES ISOFORM COUNT MATRIX CONVERTER--------------------------
geneMatrixFilePath_day55 = "/data/gpfs/projects/punim0646/manveer/lr-cortdiff-day55_gene_count.csv"
count_matrix_day55 <- read.csv(geneMatrixFilePath_day55, row.names = 1)
View(count_matrix_day55)
geneMatrixFilePath_day25 = "/data/gpfs/projects/punim0646/manveer/lr-cortdiff-day25_gene_count.csv"
count_matrix_day25 <- read.csv(geneMatrixFilePath_day25, row.names = 1)
View(count_matrix_day25)

trans_gene_dict_day55 <- count_matrix_day55[, 1, drop = FALSE]
trans_gene_dict_day25 <- count_matrix_day25[, 1, drop = FALSE]
trans_gene_dict_day25$transcript_id <- gsub("(\\d)([BE])", "\\1,\\2", trans_gene_dict_day25$transcript_id)
trans_gene_dict_day55$transcript_id <- gsub("(\\d)([BE])", "\\1,\\2", trans_gene_dict_day55$transcript_id)

# Merge the two data frames by row names
combined_trans_gene_dict <- merge(trans_gene_dict_day25, trans_gene_dict_day55, by = "row.names", all = TRUE)
rownames(combined_trans_gene_dict) <- combined_trans_gene_dict$Row.names
combined_trans_gene_dict$Row.names <- NULL

# Define a function to combine unique comma-separated transcript IDs
combine_transcripts <- function(x, y) {
  combined <- unique(unlist(strsplit(paste(x, y, sep = ","), ",")))
  return(paste(combined, collapse = ","))
}

# Apply the combine_transcripts function to each row
final_transcript_gene_dictionary <- data.frame(
  transcript_id = mapply(combine_transcripts, combined_trans_gene_dict$transcript_id.x, combined_trans_gene_dict$transcript_id.y),
  row.names = rownames(combined_trans_gene_dict)
)

# View the result, add gene names, and save as an rds object to use in other scripts
GeneNames <- ListofGeneIDstoNames(rownames(final_transcript_gene_dictionary), gtfFilePath)
final_transcript_gene_dictionary$geneName <- GeneNames
View(final_transcript_gene_dictionary)
saveRDS(final_transcript_gene_dictionary, file = "transcript_gene_dictionary.rds")

## Read in filtered seurat objects for both samples with short reads and long reads---------------
integrated.lr.filepath <- "Relabeled_LR_days25+55-script3.rds"
integrated.sr.filepath <- "Relabeled_SR_days25+55-script3.rds"
integrated.lr.seurat.obj <- readRDS(file = integrated.lr.filepath)
integrated.sr.seurat.obj <- readRDS(file = integrated.sr.filepath)

integrated.lr.ISOLVL.filepath <- "INTEGRATED_LR_days25+55_ISOLVL.rds"
integrated.lr.ISOLVL.seurat.obj <- readRDS(file = integrated.lr.ISOLVL.filepath)

## REMOVE COLS IN METADATA THAT AREN'T NEEDED ANYMORE------------------------
integrated.lr.seurat.obj@meta.data <- integrated.lr.seurat.obj@meta.data %>%
  select(-contains("RNA_snn_res"))
integrated.sr.seurat.obj@meta.data <- integrated.sr.seurat.obj@meta.data %>%
  select(-contains("RNA_snn_res"))
integrated.lr.ISOLVL.seurat.obj@meta.data <- integrated.lr.ISOLVL.seurat.obj@meta.data %>% 
  select(-contains("RNA_snn_res"))

# Add cell labels from gene level long read seurat object to its isoform level object
identical(row.names(integrated.lr.seurat.obj@meta.data), row.names(integrated.lr.ISOLVL.seurat.obj@meta.data)) # Confirm that cell barcode order is the same
integrated.lr.ISOLVL.seurat.obj$customclassif <- integrated.lr.seurat.obj$customclassif



# Create function to get conserved markers for any given cluster----------
get_conserved.LR <- function(cluster){
  cluster <- 0
  View(integrated.lr.seurat.obj@meta.data)
  FindConservedMarkers(integrated.lr.seurat.obj,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
get_cluster_top10_GeneMarkers.LR <- function(cluster, findall = FALSE){
  conserved.markers <- map_dfr(cluster, get_conserved.LR)

  top10 <- conserved.markers %>% 
    mutate(avg_fc = (LR_DAY25_avg_log2FC + LR_DAY55_avg_log2FC) /2) %>% 
    group_by(cluster_id) %>% 
    top_n(n = 10,
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
                                  top=textGrob(paste0("Gene Markers of Interest - Cluster ", cluster,"\n"),
                                               gp = gpar(fontsize = 12, fontface = "bold")))
  
  grid.layout.FoI
}
generate.VlnPlots.FOI.LR <- function(seurat.obj, FOIs, cluster){
  FoI_figures <- list()
  rows <- ceiling(length(FOIs)/2)
  
  # Convert gene IDs in Features of Interest (FOIs) list to their respective symbols using ListofGeneIDstoNames() function
  converted_geneNames <- ListofGeneIDstoNames(FOIs, gtfFilePath = gtfFilePath)
  
  for (index in seq_along(FOIs)) {
    feature_name <- converted_geneNames[index]
    FoI_figures[[index]] <- VlnPlot(object = seurat.obj, 
                                        features = FOIs[[index]]) +
      labs(title = feature_name) + NoLegend()
  }
  
  grid.layout.FoI <- grid.arrange(grobs = FoI_figures, nrow = rows, ncol = 2,
                                  top=textGrob(paste0("Gene Markers of Interest - Cluster ", cluster,"\n"),
                                               gp = gpar(fontsize = 12, fontface = "bold")))
  
  grid.layout.FoI
}

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

### Revert some of the labels to re-examine scType generated labels
#integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "4"] <- "Schwann precursor cells"
#integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "3"] <- "Microglia"
#integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "7"] <- "Immature neurons"
#integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "0"] <- "Mature neurons"
#integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "9"] <- "Endothelial cells"
#integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "6"] <- "Radial glial cells"

#integrated.sr.seurat.obj$customclassif[integrated.sr.seurat.obj$seurat_clusters == "6"] <- "Unknown"
#integrated.sr.seurat.obj$customclassif[integrated.sr.seurat.obj$seurat_clusters == "3"] <- "Mature neurons"

## INVESTIGATE CLUSTER 8 - DGE--------------------------------------
# Extract top 10 markers per cluster
cluster8_top10 <- get_cluster_top10_GeneMarkers.LR(8)

View(cluster8_top10)

# Plot interesting marker gene expression for cluster 8 (SLC17A6, STMN2, STMN4, GAP43, INA)
cluster8_FoIs <- list("ENSG00000091664.9", "ENSG00000104435.14", 
                   "ENSG00000015592.17", "ENSG00000172020.13", 
                   "ENSG00000148798.11")

feature.plts.FoIs.cluster8 <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj, cluster8_FoIs, 8)

FeaturePlot(object = integrated.lr.seurat.obj, 
            features = cluster8_FoIs,
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE,
            split.by = 'orig.ident')

vln.plts.FoIs.cluster8 <- generate.VlnPlots.FOI.LR(integrated.lr.seurat.obj, cluster8_FoIs, 8)

feature.plts.FoIs.cluster8
## INVESTIGATE CLUSTERS 9 AND 6 (hypothesized to be neural stem cells)-----------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$seurat_clusters
cluster9_top10 <- get_cluster_top10_GeneMarkers.LR(9)
View(cluster9_top10)
clusters6_top10 <- get_cluster_top10_GeneMarkers.LR(6)
View(clusters6_top10)


lineage.markers <- ListofGeneNamestoIDs(list("EOMES", "TBR1", "EMX2", "GLI3", "ID4", "HES1"), gtfFilePath)
FeaturePlot(integrated.lr.seurat.obj, reduction = 'umap', features = lineage.markers)
lineage.markers

Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$customclassif
clusters_NeuralStemCells <- get_cluster_top10_GeneMarkers.LR('Neural stem cells')
View(clusters_NeuralStemCells)

neuralstemcells_FoIs <- list("ENSG00000162493.17", "ENSG00000198467.16", "ENSG00000121966.8",
                             "ENSG00000118523.6", "ENSG00000116729.14")
feature.plts.FoIs.9AND6.cluster <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj, neuralstemcells_FoIs,
                                                            "9 + 6")
vln.plts.FoIs.9AND6.cluster <- generate.VlnPlots.FOI.LR(integrated.lr.seurat.obj, neuralstemcells_FoIs, "9 + 6")
grid.draw(feature.plts.FoIs.9AND6.cluster)

# There is more evidence supporting the notion that these are a subtype of neural progenitor cells
integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "6"] <- "Neural Progenitor cells"
integrated.lr.seurat.obj$customclassif[integrated.lr.seurat.obj$seurat_clusters == "9"] <- "Neural Progenitor cells"

## WHATS DIFFERENTIALLY EXPRESSED BETWEEN CLUSTER 8 AND CLUSTER 5
## INVESTIGATE CLUSTER 4 + 3 + 1 (hypothesized to be radial glia)--------------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$seurat_clusters
cluster4_top10 <- get_cluster_top10_GeneMarkers.LR(4)
View(cluster4_top10)
cluster3_top10 <- get_cluster_top10_GeneMarkers.LR(3)
View(clusters3_top10)
cluster1_top10 <- get_cluster_top10_GeneMarkers.LR(1)
View(cluster1_top10)

Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$customclassif

radialgliacluster_top10 <- get_cluster_top10_GeneMarkers.LR("Radial glial cells")
View(radialgliacluster_top10)

radialgliacluster_FoIs <- list("ENSG00000172201.12", "ENSG00000007372.25", 
                               "ENSG00000170370.12", "ENSG00000286522.2",
                               "ENSG00000197409.8")
feature.plts.FoIs.clusters.4.1.3 <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj,
                                                       radialgliacluster_FoIs, "4 + 3 + 1")
vln.plts.FoIs.clusters.4.1.3 <- generate.VlnPlots.FOI.LR(integrated.lr.seurat.obj,
                                                         radialgliacluster_FoIs, "4 + 3 + 1")

feature.plts.FoIs.cluster4 <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj,
                                                       cluster4_top10$gene, "4")
feature.plts.FoIs.cluster3 <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj, 
                                                       clusters3_top10$gene, "3")
feature.plts.FoIs.cluster1 <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj, 
                                                      cluster1_top10$gene, "1")

## INVESTIGATE CLUSTER 7 (hypothesized to be glutamatergic neurons)-------------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$seurat_clusters
glutamatergic_top10 <- get_cluster_top10_GeneMarkers.LR("Glutamatergic neurons")
View(cluster7_top10)

cluster7_FoIs <- list("ENSG00000112186.13", "ENSG00000105855.10",
                      "ENSG00000105088.9", "ENSG00000082684.16")
feature.plts.FoIs.cluster7 <- generate.DimPlots.FOI.LR(integrated.lr.seurat.obj,
                                                       cluster7_FoIs, "7")

View(clusters0_top10)
### Investigate markers for cluster 0, 7, and 8--------------------------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$seurat_clusters
View(head(LR.Markers_C0, 10))
top5_avg.log2FC_C0 <- head(LR.Markers_C0, 10)
FeaturePlot(integrated.lr.seurat.obj, features = top5_avg.log2FC_C0$gene,
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
top5_avg.log2FC_C0$geneName


View(head(LR.Markers_C7, 10))
top5_avg.log2FC_C7 <- head(LR.Markers_C7, 10)
FeaturePlot(integrated.lr.seurat.obj, features = top5_avg.log2FC_C7$gene,
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)


top5_avg.log2FC_C8 <- head(LR.Markers_C8, 10)
View(head(LR.Markers_C8, 10))
top5_avg.log2FC_C8 <- head(LR.Markers_C8, 10)
FeaturePlot(integrated.lr.seurat.obj, features = top5_avg.log2FC_C8$gene,
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

# Find differences between cluster 8 versus clusters 7 and 0 which I think should be merged-----------------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$customclassif
C7andC0.vs.C8 <- FindMarkers(integrated.lr.seurat.obj, ident.1 = 'Mature neurons', ident.2 = 'Glutamatergic neurons')
C7andC0.vs.C8$geneName <- ListofGeneIDstoNames(row.names(C7andC0.vs.C8), gtfFilePath)

C7andC0.vs.C8 <- C7andC0.vs.C8  %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

top15_avg.log2FC_C7andC0.vs.C8 <- head(C7andC0.vs.C8, 15)
View(top15_avg.log2FC_C7andC0.vs.C8)
FeaturePlot(integrated.lr.seurat.obj, features = row.names(top15_avg.log2FC_C7andC0.vs.C8),
            order = TRUE,
            min.cutoff = 'q10')



GABAergic.vs.C8 <- FindMarkers(integrated.lr.seurat.obj, ident.1 = 'Mature neurons', ident.2 = 'GABAergic neurons')
GABAergic.vs.C8$geneName <- ListofGeneIDstoNames(row.names(GABAergic.vs.C8), gtfFilePath)

GABAergic.vs.C8 <- GABAergic.vs.C8  %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

top15_avg.log2FC_GABAergic.vs.C8 <- head(GABAergic.vs.C8, 15)
View(top15_avg.log2FC_GABAergic.vs.C8)
FeaturePlot(integrated.lr.seurat.obj, features = row.names(top15_avg.log2FC_GABAergic.vs.C8),
            order = TRUE,
            min.cutoff = 'q10')




### Custom heatmap function using Complex Heatmap for de novo clustering -------------------------
DoClusteredHeatmap <- function(seurat.obj, markersOfInterest, longRead = FALSE, 
                               longRead.IsoLvl = FALSE, bySeuratCluster = TRUE){
  ## Some preliminary data wrangling to get a matrix for Complex Heatmap to work with
  if(bySeuratCluster){
    Idents(seurat.obj) <- seurat.obj$seurat_clusters    
  } else{
    Idents(seurat.obj) <- seurat.obj$customclassif
  }
  
  seurat.obj.counts.log <- as.data.frame(GetAssayData(seurat.obj, slot = 'scale.data'))
  if(bySeuratCluster){
    cell_clusters <- as.data.frame(seurat.obj$seurat_clusters)
    cell_order <- cell_clusters %>%
      arrange('seurat.obj$seurat_clusters')    
  }else{
    cell_clusters <- as.data.frame(seurat.obj$customclassif)
    cell_order <- cell_clusters %>% 
      arrange('seurat.obj$customclassif')
    cell_order$`seurat.obj$customclassif` <-  gsub(" ", "\n", cell_order$`seurat.obj$customclassif`)
  }

  # get rid of repeated genes in our combined top10Markers.LR list
  markersOfInterest <- markersOfInterest %>% 
    distinct(gene, .keep_all = TRUE)
  
  counts.log.top10 <- seurat.obj.counts.log %>% 
    # Filter top 10 gene markers
    rownames_to_column("gene") %>% 
    filter(gene %in% markersOfInterest$gene) %>% 
    # Remove outliers
    mutate_if(is.numeric, ~ifelse(. > 2.5, 2.5, .)) %>%
    #Reorder cells
    select(gene, all_of(rownames(cell_order))) %>% 
    # Reorder genes
    mutate(gene = factor(gene, levels = markersOfInterest$gene)) %>% 
    arrange(gene) %>% 
    # Convert to matrix (complex heatmap only uses matrices)
    column_to_rownames('gene') %>% 
    as.matrix
  
  ## Moving on to actually building the complex heatmap
  # Make colour gradient to span the feature expression values
  # Force the colour gradient to center at 0.5 to match Seurat
  #col_min <- min(counts.log.top10, na.rm = TRUE)
  #col_max <- max(counts.log.top10, na.rm = TRUE)
  col_min <- -2
  col_max <- 2
  # Create colour gradient
  custom_blue <- rgb(58, 112, 157, maxColorValue = 255)
  custom_red <- rgb(174, 40, 50, maxColorValue = 255)
  full_col_gradient <- colorRamp2(c(col_min, 0, col_max),
                                  c(custom_blue, "white", custom_red))
  
  # format cell cluster annotation (at the top of the heatmap)
  if(bySeuratCluster){
    heatmap_clusters <- cell_order$'seurat.obj$seurat_clusters'   
  } else{
    heatmap_clusters <- cell_order$'seurat.obj$customclassif'
  }
  names(heatmap_clusters) <- rownames(cell_order)  

  # Make the cluster colours the same as the ggplot default (used by seurat)
  heatmap_colours <- scales::hue_pal()(10)
  names(heatmap_colours) <- as.character(c(0:9))
  if(!bySeuratCluster){
    names(heatmap_colours) <- as.character(unique(cell_order$`seurat.obj$customclassif`))
    heatmap_colours <- heatmap_colours[as.character(unique(cell_order$`seurat.obj$customclassif`))]
  }

  # Format ComplexHeatmap top annotation
  cluster_annotation <- HeatmapAnnotation(
    Identity = heatmap_clusters,
    col = list(Identity = heatmap_colours),
    show_annotation_name = FALSE)
  
  # Convert rownames from geneids to gene names if dealing with long reads
  View(counts.log.top10)
  if(longRead){
    row.names(counts.log.top10) <- ListofGeneIDstoNames(row.names(counts.log.top10),
                                                        gtfFilePath)
  }
  if(longRead.IsoLvl){ 
    row.names(markersOfInterest) <- markersOfInterest$gene
    
    # Sort row names of our count matrix by alphabetical order
    counts.log.top10 <- counts.log.top10[order(rownames(counts.log.top10)),]
    markersOfInterest <- markersOfInterest[order(rownames(markersOfInterest)),]
    row.names(counts.log.top10) <- markersOfInterest$geneName
    View(counts.log.top10)
  }
  # Make heatmap
  Heatmap_de_novo_clustOrder <- Heatmap(counts.log.top10, use_raster = FALSE,
                                        cluster_rows = TRUE, cluster_columns = TRUE,
                                        name = "Expression \n Level",
                                        column_split = heatmap_clusters,
                                        col = full_col_gradient,
                                        show_column_names = FALSE,
                                        row_names_side = "left",
                                        top_annotation = cluster_annotation,
                                        row_names_gp = gpar(fontsize = 6))
  
  Heatmap_de_novo_clustOrder
}

### FIND ALL MARKERS LONG READ DATA (UNSUPERVISED CLUSTERS)--------------------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$seurat_clusters

lr.all.markers.only.pos <- FindAllMarkers(integrated.lr.seurat.obj,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
lr.all.markers.pos.AND.neg <- FindAllMarkers(integrated.lr.seurat.obj,
                                 logfc.threshold = 0.25,
                                 min.pct = 0.1)

#lr.all.markers <- lr.all.markers.pos.AND.neg
lr.all.markers <- lr.all.markers.only.pos

lr.all.markers$geneName <- ListofGeneIDstoNames(lr.all.markers$gene, gtfFilePath)

num.Genes <- 10

LR.Markers_C0 <- lr.all.markers %>%
  dplyr::filter(cluster == 0 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C0 <- head(LR.Markers_C0, num.Genes)

LR.Markers_C1 <- lr.all.markers %>%
  dplyr::filter(cluster == 1 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C1 <- head(LR.Markers_C1, num.Genes)

LR.Markers_C2 <- lr.all.markers %>%
  dplyr::filter(cluster == 2 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C2 <- head(LR.Markers_C2, num.Genes)

LR.Markers_C3 <- lr.all.markers %>%
  dplyr::filter(cluster == 3 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C3 <- head(LR.Markers_C3, num.Genes)

LR.Markers_C4 <- lr.all.markers %>%
  dplyr::filter(cluster == 4 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C4 <- head(LR.Markers_C4, num.Genes)

LR.Markers_C5 <- lr.all.markers %>%
  dplyr::filter(cluster == 5 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C5 <- head(LR.Markers_C5, num.Genes)

LR.Markers_C6 <- lr.all.markers %>%
  dplyr::filter(cluster == 6 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C6 <- head(LR.Markers_C6, num.Genes)

LR.Markers_C7 <- lr.all.markers %>%
  dplyr::filter(cluster == 7 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C7 <- head(LR.Markers_C7, num.Genes)

LR.Markers_C8 <- lr.all.markers %>%
  dplyr::filter(cluster == 8 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C8 <- head(LR.Markers_C8, num.Genes)

#LR.Markers_C9 <- lr.all.markers %>%
#  dplyr::filter(cluster == 9 & p_val_adj < 0.05) %>%
#  arrange(desc(abs(avg_log2FC)))
#top.n.LR.Markers_C9 <- head(LR.Markers_C9, num.Genes)


combined.top10Markers.LR <- bind_rows(top.n.LR.Markers_C0, top.n.LR.Markers_C1, top.n.LR.Markers_C2,
                                      top.n.LR.Markers_C3, top.n.LR.Markers_C4, top.n.LR.Markers_C5,
                                      top.n.LR.Markers_C6, top.n.LR.Markers_C7, top.n.LR.Markers_C8)

### FIND ALL MARKERS SHORT READ DATA (UNSUPERVISED CLUSTERS)--------------------------
Idents(integrated.sr.seurat.obj) <- integrated.sr.seurat.obj$seurat_clusters
sr.all.markers.only.pos <- FindAllMarkers(integrated.sr.seurat.obj,
                                          logfc.threshold = 0.25,
                                          min.pct = 0.1,
                                          only.pos = TRUE)
sr.all.markers.pos.AND.neg <- FindAllMarkers(integrated.sr.seurat.obj,
                                             logfc.threshold = 0.25,
                                             min.pct = 0.1)

#SR.all.markers <- sr.all.markers.pos.AND.neg
SR.all.markers <- sr.all.markers.only.pos


num.Genes <- 10

SR.Markers_C0 <- SR.all.markers %>%
  dplyr::filter(cluster == 0 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C0 <- head(SR.Markers_C0, num.Genes)

SR.Markers_C1 <- SR.all.markers %>%
  dplyr::filter(cluster == 1 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C1 <- head(SR.Markers_C1, num.Genes)

SR.Markers_C2 <- SR.all.markers %>%
  dplyr::filter(cluster == 2 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C2 <- head(SR.Markers_C2, num.Genes)

SR.Markers_C3 <- SR.all.markers %>%
  dplyr::filter(cluster == 3 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C3 <- head(SR.Markers_C3, num.Genes)

SR.Markers_C4 <- SR.all.markers %>%
  dplyr::filter(cluster == 4 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C4 <- head(SR.Markers_C4, num.Genes)

SR.Markers_C5 <- SR.all.markers %>%
  dplyr::filter(cluster == 5 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C5 <- head(SR.Markers_C5, num.Genes)

SR.Markers_C6 <- SR.all.markers %>%
  dplyr::filter(cluster == 6 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C6 <- head(SR.Markers_C6, num.Genes)

SR.Markers_C7 <- SR.all.markers %>%
  dplyr::filter(cluster == 7 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C7 <- head(SR.Markers_C7, num.Genes)

SR.Markers_C8 <- SR.all.markers %>%
  dplyr::filter(cluster == 8 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C8 <- head(SR.Markers_C8, num.Genes)

SR.Markers_C9 <- SR.all.markers %>%
  dplyr::filter(cluster == 9 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.SR.Markers_C9 <- head(SR.Markers_C9, num.Genes)

combined.top10Markers.SR <- bind_rows(top.n.SR.Markers_C0, top.n.SR.Markers_C1, top.n.SR.Markers_C2,
                                      top.n.SR.Markers_C3, top.n.SR.Markers_C4, top.n.SR.Markers_C5,
                                      top.n.SR.Markers_C6, top.n.SR.Markers_C7, top.n.SR.Markers_C8,
                                      top.n.SR.Markers_C9)




### FIND ALL MARKERS LONG READ DATA (LABELED CELL TYPES)---------------------------------
Idents(integrated.lr.seurat.obj) <- integrated.lr.seurat.obj$customclassif
CellType.lr.all.markers.only.pos <- FindAllMarkers(integrated.lr.seurat.obj,
                                          logfc.threshold = 0.25,
                                          min.pct = 0.1,
                                          only.pos = TRUE)
CellType.lr.all.markers.pos.AND.neg <- FindAllMarkers(integrated.lr.seurat.obj,
                                             logfc.threshold = 0.25,
                                             min.pct = 0.1)

#CellType.lr.all.markers <- CellType.lr.all.markers.pos.AND.neg
CellType.lr.all.markers <- CellType.lr.all.markers.only.pos

CellType.lr.all.markers$geneName <- ListofGeneIDstoNames(CellType.lr.all.markers$gene, gtfFilePath)

LR.cellType.labels <- c("Neural Progenitor cells", "Mature neurons",
                     "Radial glial cells", "Endothelial cells",
                     "Myelinating Schwann cells", "Immature neurons")

num.Genes <- 10

CellType.LR.Markers_C1 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[1] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C1 <- head(CellType.LR.Markers_C1, num.Genes)

CellType.LR.Markers_C2 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[2] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C2 <- head(CellType.LR.Markers_C2, num.Genes)

CellType.LR.Markers_C3 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[3] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C3 <- head(CellType.LR.Markers_C3, num.Genes)

CellType.LR.Markers_C4 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[4] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C4 <- head(CellType.LR.Markers_C4, num.Genes)

CellType.LR.Markers_C5 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[5] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C5 <- head(CellType.LR.Markers_C5, num.Genes)

CellType.LR.Markers_C6 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[6] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C6 <- head(CellType.LR.Markers_C6, num.Genes)

CellType.LR.Markers_C7 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[7] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C7 <- head(CellType.LR.Markers_C7, num.Genes)

CellType.LR.Markers_C8 <- CellType.lr.all.markers %>%
  dplyr::filter(cluster == LR.cellType.labels[8] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C8 <- head(CellType.LR.Markers_C8, num.Genes)


combined.top10Markers.LR.CellType <- bind_rows(CellType.top.n.LR.Markers_C1, CellType.top.n.LR.Markers_C2, 
                                      CellType.top.n.LR.Markers_C3, CellType.top.n.LR.Markers_C4,
                                      CellType.top.n.LR.Markers_C5, CellType.top.n.LR.Markers_C6)


### FIND ALL MARKERS SHORT READ DATA (LABELED CELL TYPES)---------------------------------
Idents(integrated.sr.seurat.obj) <- integrated.sr.seurat.obj$customclassif
CellType.SR.all.markers.only.pos <- FindAllMarkers(integrated.sr.seurat.obj,
                                                   logfc.threshold = 0.25,
                                                   min.pct = 0.1,
                                                   only.pos = TRUE)
CellType.SR.all.markers.pos.AND.neg <- FindAllMarkers(integrated.sr.seurat.obj,
                                                      logfc.threshold = 0.25,
                                                      min.pct = 0.1)

#CellType.SR.all.markers <- CellType.SR.all.markers.pos.AND.neg
CellType.SR.all.markers <- CellType.SR.all.markers.only.pos

SR.cellType.labels <- c("Cancer cells", "GABAergic neurons",
                        "Endothelial cells", "Radial glial cells",
                        "Neuroblasts", "Mature neurons",
                        "Neural Progenitor cells", "Immature neurons")

num.Genes <- 10

CellType.SR.Markers_C1 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[1] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C1 <- head(CellType.SR.Markers_C1, num.Genes)

CellType.SR.Markers_C2 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[2] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C2 <- head(CellType.SR.Markers_C2, num.Genes)

CellType.SR.Markers_C3 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[3] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C3 <- head(CellType.SR.Markers_C3, num.Genes)

CellType.SR.Markers_C4 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[4] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C4 <- head(CellType.SR.Markers_C4, num.Genes)

CellType.SR.Markers_C5 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[5] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C5 <- head(CellType.SR.Markers_C5, num.Genes)

CellType.SR.Markers_C6 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[6] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C6 <- head(CellType.SR.Markers_C6, num.Genes)

CellType.SR.Markers_C7 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[7] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C7 <- head(CellType.SR.Markers_C7, num.Genes)

CellType.SR.Markers_C8 <- CellType.SR.all.markers %>%
  dplyr::filter(cluster == SR.cellType.labels[8] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.SR.Markers_C8 <- head(CellType.SR.Markers_C8, num.Genes)

combined.top10Markers.SR.CellType <- bind_rows(CellType.top.n.SR.Markers_C1, CellType.top.n.SR.Markers_C2, 
                                               CellType.top.n.SR.Markers_C3, CellType.top.n.SR.Markers_C4,
                                               CellType.top.n.SR.Markers_C5, CellType.top.n.SR.Markers_C6,
                                               CellType.top.n.SR.Markers_C7, CellType.top.n.SR.Markers_C8)





### FIND ALL MARKERS LONG READ DATA (ISOFORM LEVEL) (UNSUPERVISED CLUSTERS)--------------------------
Idents(integrated.lr.ISOLVL.seurat.obj) <- integrated.lr.ISOLVL.seurat.obj$seurat_clusters

lr.all.markers.only.pos.ISOLVL <- FindAllMarkers(integrated.lr.ISOLVL.seurat.obj,
                                          logfc.threshold = 0.25,
                                          min.pct = 0.1,
                                          only.pos = TRUE)
lr.all.markers.pos.AND.neg.ISOLVL <- FindAllMarkers(integrated.lr.ISOLVL.seurat.obj,
                                             logfc.threshold = 0.25,
                                             min.pct = 0.1)

#lr.all.markers <- lr.all.markers.pos.AND.neg
lr.all.markers.ISOLVL <- lr.all.markers.only.pos.ISOLVL
View(lr.all.markers.ISOLVL)

lr.all.markers.ISOLVL$geneName <- ListofTranscriptIDstoNames(lr.all.markers.ISOLVL$gene, gtfFilePath)

num.Genes <- 10

LR.Markers_C0.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 0 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C0.ISOLVL <- head(LR.Markers_C0.ISOLVL, num.Genes)

LR.Markers_C1.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 1 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C1.ISOLVL <- head(LR.Markers_C1.ISOLVL, num.Genes)

LR.Markers_C2.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 2 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C2.ISOLVL <- head(LR.Markers_C2.ISOLVL, num.Genes)

LR.Markers_C3.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 3 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C3.ISOLVL <- head(LR.Markers_C3.ISOLVL, num.Genes)

LR.Markers_C4.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 4 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C4.ISOLVL <- head(LR.Markers_C4.ISOLVL, num.Genes)

LR.Markers_C5.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 5 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C5.ISOLVL <- head(LR.Markers_C5.ISOLVL, num.Genes)

LR.Markers_C6.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 6 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C6.ISOLVL <- head(LR.Markers_C6.ISOLVL, num.Genes)

LR.Markers_C7.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 7 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C7.ISOLVL <- head(LR.Markers_C7.ISOLVL, num.Genes)

LR.Markers_C8.ISOLVL <- lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == 8 & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
top.n.LR.Markers_C8.ISOLVL <- head(LR.Markers_C8.ISOLVL, num.Genes)

#LR.Markers_C9 <- lr.all.markers %>%
#  dplyr::filter(cluster == 9 & p_val_adj < 0.05) %>%
#  arrange(desc(abs(avg_log2FC)))
#top.n.LR.Markers_C9 <- head(LR.Markers_C9, num.Genes)


combined.top10Markers.LR.ISOLVL <- bind_rows(top.n.LR.Markers_C0.ISOLVL, top.n.LR.Markers_C1.ISOLVL, top.n.LR.Markers_C2.ISOLVL,
                                      top.n.LR.Markers_C3.ISOLVL, top.n.LR.Markers_C4.ISOLVL, top.n.LR.Markers_C5.ISOLVL,
                                      top.n.LR.Markers_C6.ISOLVL, top.n.LR.Markers_C7.ISOLVL, top.n.LR.Markers_C8.ISOLVL)

View(combined.top10Markers.LR.ISOLVL)

combined.top10Markers.LR.ISOLVL[8]

## Retrieve the genes for each novel Bambu isoform, add that information into the geneName column
for (i in seq_along(combined.top10Markers.LR.ISOLVL$gene)) {
  gene <- combined.top10Markers.LR.ISOLVL$gene[i]
  if (grepl("Bambu", gene)) {
    matches <- grepl(gene, final_transcript_gene_dictionary$transcript_id)
    matched_gene_id <- rownames(final_transcript_gene_dictionary)[matches]
    
    # Assuming there's only one match
    convertedGeneName <- FeatureIDtoName(matched_gene_id[1], ConvType = "g.to.gName", gtfFilePath = gtfFilePath)
    combined.top10Markers.LR.ISOLVL$geneName[i] <- paste0("(Novel ", convertedGeneName, " Isoform) ", gene)
  }
}


### FIND ALL MARKERS LONG READ DATA (ISOFORM LEVEL) (LABELED CELL TYPES)---------------------------------
Idents(integrated.lr.ISOLVL.seurat.obj) <- integrated.lr.ISOLVL.seurat.obj$customclassif
CellType.lr.all.markers.only.pos.ISOLVL <- FindAllMarkers(integrated.lr.ISOLVL.seurat.obj,
                                                   logfc.threshold = 0.25,
                                                   min.pct = 0.1,
                                                   only.pos = TRUE)
CellType.lr.all.markers.pos.AND.neg.ISOLVL <- FindAllMarkers(integrated.lr.ISOLVL.seurat.obj,
                                                      logfc.threshold = 0.25,
                                                      min.pct = 0.1)



#CellType.lr.all.markers.ISOLVL <- CellType.lr.all.markers.pos.AND.neg.ISOLVL
CellType.lr.all.markers.ISOLVL <- CellType.lr.all.markers.only.pos.ISOLVL

CellType.lr.all.markers.ISOLVL$geneName <- ListofTranscriptIDstoNames(CellType.lr.all.markers.ISOLVL$gene, gtfFilePath)


LR.cellType.labels <- c("Neural Progenitor cells", "Mature neurons",
                        "Radial glial cells", "Endothelial cells",
                        "Myelinating Schwann cells", "Immature neurons")

num.Genes <- 10

CellType.LR.Markers_C1.ISOLVL <- CellType.lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == LR.cellType.labels[1] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C1.ISOLVL <- head(CellType.LR.Markers_C1.ISOLVL, num.Genes)

CellType.LR.Markers_C2.ISOLVL <- CellType.lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == LR.cellType.labels[2] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C2.ISOLVL <- head(CellType.LR.Markers_C2.ISOLVL, num.Genes)

CellType.LR.Markers_C3.ISOLVL <- CellType.lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == LR.cellType.labels[3] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C3.ISOLVL <- head(CellType.LR.Markers_C3.ISOLVL, num.Genes)

CellType.LR.Markers_C4.ISOLVL <- CellType.lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == LR.cellType.labels[4] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C4.ISOLVL <- head(CellType.LR.Markers_C4.ISOLVL, num.Genes)

CellType.LR.Markers_C5.ISOLVL <- CellType.lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == LR.cellType.labels[5] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C5.ISOLVL <- head(CellType.LR.Markers_C5.ISOLVL, num.Genes)

CellType.LR.Markers_C6.ISOLVL <- CellType.lr.all.markers.ISOLVL %>%
  dplyr::filter(cluster == LR.cellType.labels[6] & p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))
CellType.top.n.LR.Markers_C6.ISOLVL <- head(CellType.LR.Markers_C6.ISOLVL, num.Genes)



combined.top10Markers.LR.CellType.ISOLVL <- bind_rows(CellType.top.n.LR.Markers_C1.ISOLVL, CellType.top.n.LR.Markers_C2.ISOLVL, 
                                               CellType.top.n.LR.Markers_C3.ISOLVL, CellType.top.n.LR.Markers_C4.ISOLVL,
                                               CellType.top.n.LR.Markers_C5.ISOLVL, CellType.top.n.LR.Markers_C6.ISOLVL)


## Retrieve the genes for each novel Bambu isoform, add that information into the geneName column
for (i in seq_along(combined.top10Markers.LR.CellType.ISOLVL$gene)) {
  gene <- combined.top10Markers.LR.CellType.ISOLVL$gene[i]
  if (grepl("Bambu", gene)) {
    matches <- grepl(gene, final_transcript_gene_dictionary$transcript_id)
    matched_gene_id <- rownames(final_transcript_gene_dictionary)[matches]
    
    # Assuming there's only one match
    convertedGeneName <- FeatureIDtoName(matched_gene_id[1], ConvType = "g.to.gName", gtfFilePath = gtfFilePath)
    combined.top10Markers.LR.CellType.ISOLVL$geneName[i] <- paste0("(Novel ", convertedGeneName, " Isoform) ", gene)
  }
}


#### GENERATE DGE HEATMAPS-----------------------------------
lr.unsupervised.cluster.heatmap <- DoClusteredHeatmap(integrated.lr.seurat.obj, combined.top10Markers.LR, 
                                                      longRead = TRUE)
lr.cellType.cluster.heatmap <- DoClusteredHeatmap(integrated.lr.seurat.obj, combined.top10Markers.LR, 
                                                  longRead = TRUE,
                                                  bySeuratCluster = FALSE)

lr.ISOLVL.unsupervised.cluster.heatmap <- DoClusteredHeatmap(integrated.lr.ISOLVL.seurat.obj, combined.top10Markers.LR.ISOLVL,
                                                          longRead.IsoLvl = TRUE)
lr.ISOLVL.celltype.cluster.heatmap <- DoClusteredHeatmap(integrated.lr.ISOLVL.seurat.obj, 
                                                         combined.top10Markers.LR.CellType.ISOLVL,
                                                         longRead.IsoLvl = TRUE,
                                                         bySeuratCluster = FALSE)


sr.unsupervised.cluster.heatmap <- DoClusteredHeatmap(integrated.sr.seurat.obj, combined.top10Markers.SR, 
                                                      longRead = FALSE)
sr.cellType.cluster.heatmap <- DoClusteredHeatmap(integrated.sr.seurat.obj, combined.top10Markers.SR,
                                                  longRead = FALSE,
                                                  bySeuratCluster = FALSE)

sr.cellType.cluster.heatmap.wCorrectMarkers <- DoClusteredHeatmap(integrated.sr.seurat.obj,
                                                                 combined.top10Markers.SR.CellType,
                                                                 longRead = FALSE,
                                                                 bySeuratCluster = FALSE)
lr.cellType.cluster.heatmap.wCorrectMarkers <- DoClusteredHeatmap(integrated.lr.seurat.obj,
                                                                combined.top10Markers.LR.CellType,
                                                                longRead = TRUE,
                                                                bySeuratCluster = FALSE)


lr.unsupervised.cluster.heatmap
lr.cellType.cluster.heatmap.wCorrectMarkers

lr.ISOLVL.unsupervised.cluster.heatmap
lr.ISOLVL.celltype.cluster.heatmap

sr.unsupervised.cluster.heatmap
sr.cellType.cluster.heatmap.wCorrectMarkers


### Export our DGE heatmaps derived from FindAllMarkers() function----------------------
pdf("LongRead-heatmaps.pdf", width = 11, height = 8)
grid.draw(lr.celltype | lr.unsupervised)
draw(lr.unsupervised.cluster.heatmap)
draw(lr.ISOLVL.unsupervised.cluster.heatmap)

draw(lr.cellType.cluster.heatmap.wCorrectMarkers)
draw(lr.ISOLVL.celltype.cluster.heatmap)
dev.off()


pdf("ShortRead-heatmaps.pdf", width = 11, height = 8)
grid.draw(sr.celltype | sr.unsupervised)
draw(sr.unsupervised.cluster.heatmap)

draw(sr.cellType.cluster.heatmap.wCorrectMarkers)
dev.off()


### Prepare files for GO enrichment analysis------------------------------------------------------
exportBackgroundGeneSet <- function(seurat.obj, longRead = FALSE, csvFileName){
  gene_list <- rownames(seurat.obj@assays$RNA@counts)
  gene_df <- data.frame(Gene = gene_list)
  
  if(longRead){
    # Remove rows where Gene column contains "Bambu"
    gene_df <- gene_df[!grepl("Bambu", gene_df$Gene), , drop = FALSE]
    # Remove substring "-PAR-Y" from the Gene column
    gene_df$Gene <- gsub("-PAR-Y", "", gene_df$Gene)
    # Convert Gene IDs to their corresponding Gene Names
    gene_df <- as.data.frame(ListofGeneIDstoNames(gene_df$Gene, gtfFilePath))
  }
  
  # Write to CSV
  write.csv(gene_df, csvFileName, row.names = FALSE)
}

# Create csv files that contain the top 10 genes for GO enrichment analysis
# (sorted by avg log fold change, after filtering for genes with an adjusted p value <= 0.05)
nested.list.LR <- split(lr.all.markers, lr.all.markers$cluster)
nested.list.SR <- split(SR.all.markers, SR.all.markers$cluster)
nested.list.LR.CellType <- split(CellType.lr.all.markers, 
                               CellType.lr.all.markers$cluster)
nested.list.SR.CellType <- split(CellType.SR.all.markers,
                                 CellType.SR.all.markers$cluster)

exportDEGsAsCSV <- function(nested.list, sequencingType, longRead = FALSE){
  for(i in 1:length(nested.list)){
    csvFileName <- paste0(sequencingType,"_Cluster_", nested.list[[i]]$cluster[1], ".csv")
    if(longRead){
      gene_df <- nested.list[[i]]$geneName
      gene_df <- gene_df[sapply(gene_df, length) > 0]
      gene_df <- as.data.frame(gene_df)
    }else{
      gene_df <- as.data.frame(nested.list[[i]]$gene)
    }
    # Write to CSV
    write.csv(gene_df, csvFileName, row.names = FALSE)
  }
}
setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline/ClusterCSVFiles-GOanalysis")
exportDEGsAsCSV(nested.list.SR, "SR")
exportDEGsAsCSV(nested.list.SR.CellType, "SR-CellType")
exportDEGsAsCSV(nested.list.LR, "LR", longRead = TRUE)
exportDEGsAsCSV(nested.list.LR.CellType, "LR-CellType", longRead = TRUE)
setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

# Tidy up csv files to remove unwanted characters : 
# Specify the directory containing the CSV files
directory_path <- "/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline/ClusterCSVFiles-GOanalysis"  # Replace with the actual path

# List all CSV files in the directory
csv_files <- list.files(path = directory_path, pattern = "\\.csv$")
# Loop through each CSV file
for (file in csv_files) {
  # Full path to the file
  full_path <- file.path(directory_path, file)
  # Read the CSV into a data frame
  df <- read.csv(full_path, stringsAsFactors = FALSE)
  # Remove unwanted characters and substrings from all columns
  df[] <- lapply(df, function(x) gsub('"', '', x))
  # Overwrite the original CSV file
  write.csv(df, full_path, row.names = FALSE, quote = FALSE)
}
# List all CSV files in the directory that start with "LR"
csv_files <- list.files(directory_path, pattern = "^LR.*\\.csv$")
# Loop through each CSV file
for (file in csv_files) {
  # Read the file, skipping the first line
  data <- read_lines(file, skip = 1)
  # Combine all lines into a single string separated by commas
  all_data <- paste(data, collapse = ",")
  # Replace all commas with newline characters
  all_data_newline <- gsub(",", "\n", all_data)
  # Write the modified content back to the file
  write_lines(all_data_newline, file)
  cat(paste("Processed", file, "\n"))
}


### GENERATE AND EXPORT FILES NEEDED FOR GO ENRICHMENT ANALYSIS------------------------
# Export our background gene sets for both long read and short read Seurat objects 
# (every single gene that is expressed in both samples)
exportBackgroundGeneSet(integrated.lr.seurat.obj, longRead = TRUE, csvFileName = "LR_background_gene_set.csv")
exportBackgroundGeneSet(integrated.sr.seurat.obj, csvFileName = "SR_background_gene_set.csv")
