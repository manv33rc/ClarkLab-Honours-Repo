# Script to perform trajectory analysis
# https://www.nature.com/articles/s41467-019=10291-0
# https://blog.bioturing.com/2022/06/13/single-cell-rna-seq-trajectory-analysis-review/
# https://bustools.github.io/BUS_notebooks_R/slingshot.html

setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

set.seed(4242)
library(reticulate)
use_python("/Users/manveerchuahan/miniconda3/bin/python3.11")
library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(slingshot)


transcript.gene.dict <- readRDS(file = 'transcript_gene_dictionary.rds')
View(transcript.gene.dict)

# Import our feature_ID_converter python script
source_python("/Volumes/Expansion/Scripts/feature_ID_converter_BACKUP.py")
gtfFilePath <- "/Volumes/Expansion/temp/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

## Read in integrated seurat objects for both samples with short reads and long reads---------------
integrated.lr.filepath <- "Relabeled_LR_days25+55-script3.rds"
integrated.sr.filepath <- "Relabeled_SR_days25+55-script3.rds"
integrated.lr.ISOLVL.filepath <- "INTEGRATED_LR_days25+55_ISOLVL.rds"

integrated.lr.ISOLVL.seurat.obj <- readRDS(file = integrated.lr.ISOLVL.filepath)
integrated.lr.seurat.obj <- readRDS(file = integrated.lr.filepath)
integrated.sr.seurat.obj <- readRDS(file = integrated.sr.filepath)

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

sr.unsupervised | sr.unsupervised.umap

DimPlot(integrated.sr.seurat.obj, reduction = 'umap', group.by = 'orig.ident')

## trajectory analysis functions (monocle3) --------------------------------------
printDEGPDF <- function(deg_list, seurat.obj, pdfFileName){
  num.pages <- ceiling(nrow(deg_list) / 9)
  nrow(deg_list)
  degs.to.print <- deg_list
  pdf.pages <- list()
  
  # Open a new PDF device
  pdf(file = pdfFileName, width = 11, height = 8)
  
  for(page in 1:num.pages){
    if(nrow(degs.to.print) >= 9){
      current.gene.queue <- rownames(degs.to.print)[1:9]
    } else{
      current.gene.queue <- rownames(degs.to.print)
    }
    pdf.pages[[page]] <- FeaturePlot(seurat.obj, features = current.gene.queue)
    
    # Print the plot to the current PDF device
    print(pdf.pages[[page]])
    
    degs.to.print <- degs.to.print[!(rownames(degs.to.print) %in% current.gene.queue), ]
    print(page)
    print(current.gene.queue)
  }
  
  # Close the PDF device
  dev.off()
  
  pdf.pages
}

# Work on seeing potential isoform switching here
performTrajectoryAnalysis <- function(seurat.obj, name, rootCluster, 
                                      byCellLabel = TRUE, generatePDF = TRUE){
  ### 1. Create a cell_data_set object------------------------------
  ## monocle3 requires a cell_data_set object, thus we need to convert our seurat object into this class
  # Seurat Wrappers has a function called as.cell_data_set() which can convert our seurat object into this class
  
  cds <- as.cell_data_set(seurat.obj)
  
  # We require a column called gene short name, lets add that column
  fData(cds)$gene_short_name <- rownames(fData(cds))
  
  
  ## 2. Cluster cells using clustering information from seurat's UMAP---------------------------------
  # assign partitions
  recreate.partition <- c(rep(1, length(cds@colData@rownames)))
  names(recreate.partition) <- cds@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  
  cds@clusters$UMAP$partitions <- recreate.partition
  
  # Assign the cluster information
  list_cluster <- seurat.obj@active.ident
  cds@clusters$UMAP$clusters <- list_cluster
  
  # Assign UMAP coordinate - cell embeddings
  cds@int_colData@listData$reducedDims$UMAP <- seurat.obj@reductions$umap@cell.embeddings
  
  # plot
  cluster.before.trajectory <- plot_cells(cds,
                                          color_cells_by = 'cluster',
                                          label_groups_by_cluster = FALSE,
                                          group_label_size = 5) +
    theme(legend.position = 'right')
  print(cluster.before.trajectory)
  
  cluster.names <- plot_cells(cds,
                              color_cells_by = 'customclassif',
                              label_groups_by_cluster = FALSE,
                              group_label_size = 5) +
    theme(legend.position = 'right') +
    scale_colour_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan', 'purple'))
  
  print(cluster.before.trajectory | cluster.names)
  
  
  ## 3. Learn trajectory graph ------------------------------
  # If use_partition = TRUE, then it will use partitions and learn a disjoint trajectory for each partition
  # since we know that all our cells follow one trajectory we set this parameter to false
  cds <- learn_graph(cds, use_partition = FALSE) 
  

  trajectory.line <- plot_cells(cds, 
             color_cells_by = 'customclassif',
             label_groups_by_cluster = FALSE,
             label_branch_points = FALSE,
             label_roots = FALSE,
             label_leaves = FALSE,
             group_label_size = 5)
  print(trajectory.line)
  
  ## 4. Order cells in pseudo-time-----------------
  # Root cells are defined by domain knowledge
  cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == rootCluster]))
  
  pseudotime.line <- plot_cells(cds,
             color_cells_by = 'pseudotime',
             label_groups_by_cluster = F,
             label_branch_points = F,
             label_roots = F,
             label_leaves = F)
  print(pseudotime.line)
  
  # Make a boxplot showing a range of pseudotimes for each cell type
  cds$monocle3_pseudotime <- pseudotime(cds)
  data.pseudo <- as.data.frame(colData(cds))
  
  if(byCellLabel){
    boxplot <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(customclassif, monocle3_pseudotime, median), fill = customclassif)) +
      geom_boxplot()
  } else{
    boxplot <- ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = seurat_clusters)) +
      geom_boxplot()
  }
  
  print(boxplot)

  
  ## 5. Find genes that change as a function of pseudo-time -----------------------------
  degs <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
  
  deg_list <- degs %>% 
    arrange(q_value) %>% 
    filter(status == 'OK') %>% 
    filter(q_value == 0)
  
  if(generatePDF){
    printDEGPDF(deg_list, seurat.obj, paste0(name,"-degs.pdf"))
  }

  ## 6. Visualizing pseudotime in seurat
  seurat.obj$pseudotime <- pseudotime(cds)
  if(byCellLabel){
    Idents(seurat.obj) <- seurat.obj$customclassif
  } else{
    Idents(seurat.obj) <- seurat.obj$seurat_clusters
  }

  pseudotime.umap <- FeaturePlot(seurat.obj, features = 'pseudotime', label = T, repel = T) + 
    labs(title = paste0(name, " Pseudotime Values"))
  
  print(pseudotime.umap)
  
  
  outputs <- list(cluster.before.trajectory, cluster.names,
                  trajectory.line, pseudotime.line,
                  boxplot, pseudotime.umap, cds)
  message("---------------TRAJECTORY ANALYSIS COMPLETE--------------------")
  
  outputs
}


## executing trajectory analysis functions---------------------------------------

# indexing outputs -> 1) cluster.before.trajectory 2) cluster.names 3) trajectory.line
# 4) pseudotime.line 5) boxplot 6) pseudotime.umap

sr.celltype | sr.unsupervised
shortread.outputs <- performTrajectoryAnalysis(seurat.obj = integrated.sr.seurat.obj, 
                                               name = "Integrated-ShortReads", 
                                               rootCluster = 8,
                                               generatePDF = FALSE)

shortread.outputs.unsupervisedClusters <- performTrajectoryAnalysis(seurat.obj = integrated.sr.seurat.obj,
                                                                    name = "Integrated-ShortReads",
                                                                    rootCluster = 8,
                                                                    byCellLabel = FALSE,
                                                                    generatePDF = FALSE)


lr.celltype | lr.unsupervised
longread.outputs <- performTrajectoryAnalysis(seurat.obj = integrated.lr.seurat.obj, 
                                              name = "Integrated-LongReads", 
                                              rootCluster = 7,
                                              generatePDF = FALSE)
lr.celltype | longread.outputs[[4]]
longread.outputs[[5]]

longread.outputs.unsupervisedClusters <- performTrajectoryAnalysis(seurat.obj = integrated.lr.seurat.obj,
                                                                   name = "Integrated-LongReads",
                                                                   rootCluster = 7,
                                                                   byCellLabel = FALSE,
                                                                   generatePDF = FALSE)
# EXPORT THESE AS PDF************ :
shortread.outputs.unsupervisedClusters[[5]] | sr.celltype | sr.unsupervised | shortread.outputs.unsupervisedClusters[[4]]
shortread.outputs[[5]] | sr.celltype | sr.unsupervised | shortread.outputs[[4]]
longread.outputs.unsupervisedClusters[[5]] | lr.celltype | lr.unsupervised | longread.outputs.unsupervisedClusters[[4]]
longread.outputs[[5]] | lr.celltype | lr.unsupervised | longread.outputs[[4]]

pdf("Short-read-Monocle3.pdf", width = 11, height = 8)
grid.draw(shortread.outputs.unsupervisedClusters[[4]])
grid.draw(sr.celltype | sr.unsupervised)
grid.draw(shortread.outputs.unsupervisedClusters[[5]] | shortread.outputs[[5]])
# Close the PDF device
dev.off()

pdf("Long-read-Monocle3.pdf", width = 11, height = 8)
grid.draw(longread.outputs.unsupervisedClusters[[4]])
grid.draw(lr.celltype | lr.unsupervised)
grid.draw(longread.outputs.unsupervisedClusters[[5]] | longread.outputs[[5]])
# Close the PDF device
dev.off()

# transfer long read gene level cell labels to its associated isoform level object
integrated.lr.ISOLVL.seurat.obj$customclassif <- integrated.lr.seurat.obj$customclassif


isolvl.output <- longread.isolevel.outputs <- performTrajectoryAnalysis(seurat.obj = integrated.lr.ISOLVL.seurat.obj, 
                                                       name = "Integrated-LongReads-ISOLVL", 
                                                       rootCluster = 7, 
                                                       byCellLabel = FALSE, 
                                                       generatePDF = TRUE)

### Trajectory analysis with slingshot------------------------------
# Code adapted from : https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html

performTrajectoryAnalysisSlingshot <- function(seurat.obj, 
                                               rootCluster, endCluster = NULL,
                                               unsupervisedUMAP, celltypeUMAP,
                                               title){
  figTitle <- paste(title, "- Slingshot Figures")
  
  # Define a color palette to use
  pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
  
  # Extract the different matrices needed for slingshot input
  dimred <- seurat.obj@reductions$umap@cell.embeddings
  clustering <- seurat.obj$seurat_clusters
  counts <- as.matrix(seurat.obj@assays$RNA@counts[seurat.obj@assays$RNA@var.features, ])
  
  # convert clustering information to the right class (from factor to numeric)
  clustering.typecast <- as.numeric(as.character(clustering))
  names(clustering.typecast) <- names(clustering)


  # Run slingshot lineage identification
  lineages <- getLineages(data = dimred, clusterLabels = clustering.typecast, 
                          start.clus = rootCluster, end.clus = endCluster)
  sds <- as.SlingshotDataSet(lineages)
  
  ## Create a pdf to store plots
  fileName <- paste0(title, "-slingshot.pdf") %>% gsub(" ", "-", .)
  
  pdf(fileName, width = 13, height = 8)
  umap.layout <- grid.arrange(unsupervisedUMAP,
                              celltypeUMAP,
                              nrow = 1, ncol = 2,
                              top=textGrob(paste0(figTitle),
                                           gp = gpar(fontsize = 12)))
  # Plot the lineages
  par(mfrow = c(1, 2))
  plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
  for (i in levels(clustering)) {
    text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
  }
  plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
  lines(sds, lwd = 3, col = "black")
  
  curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
  curves <- as.SlingshotDataSet(curves)
  plot(dimred, col = pal[clustering], asp = 1, pch = 16)
  for (i in levels(clustering)) {
    text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
  }
  lines(curves, lwd = 3, col = "black")
  dev.off()
  
  message(paste("Slingshot Trajectory Analysis completed.\n PDF has been generated :", fileName))
}

performTrajectoryAnalysisSlingshot(seurat.obj = integrated.sr.seurat.obj,
                                   rootCluster = "8",
                                   endCluster = "1",
                                   unsupervisedUMAP = sr.unsupervised,
                                   celltypeUMAP = sr.celltype,
                                   title = "Integrated Short Reads")

performTrajectoryAnalysisSlingshot(seurat.obj = integrated.lr.seurat.obj,
                                   rootCluster = "7",
                                   endCluster = "0",
                                   unsupervisedUMAP = lr.unsupervised,
                                   celltypeUMAP = lr.celltype,
                                   title = "Integrated Long Reads")

# Slingshot trajectory analysis with no end cluster defined
performTrajectoryAnalysisSlingshot(seurat.obj = integrated.lr.seurat.obj,
                                   rootCluster = "7",
                                   endCluster = NULL,
                                   unsupervisedUMAP = lr.unsupervised,
                                   celltypeUMAP = lr.celltype,
                                   title = "Integrated Long Reads No End")

performTrajectoryAnalysisSlingshot(seurat.obj = integrated.sr.seurat.obj,
                                   rootCluster = "8",
                                   endCluster = NULL,
                                   unsupervisedUMAP = sr.unsupervised,
                                   celltypeUMAP = sr.celltype,
                                   title = "Integrated Short Reads No End")
