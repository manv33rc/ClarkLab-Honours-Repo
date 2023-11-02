library(Seurat)
library(tidyverse)
library(HGNChelper)
library(igraph)
library(openxlsx)
library(ggraph)
library(data.tree)
library(stringr)
library(scales)
library(gridExtra)
library(grid)
library(reticulate)

set.seed(4242)

# load scType's gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load scType's cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# Import our feature_ID_converter python script
source_python("/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/feature_ID_converter.py")
gtfFilePath <- "/data/gpfs/projects/punim0646/manveer/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

### ASSIGN FILEPATHS TO VARIABLES------------------------------------------
sr.day25.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/SR_DAY25-QCed.rds"
sr.day55.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/SR_DAY55-QCed.rds"

q20.day25.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/Q20_FIXED-QCed.rds"
lr.day25.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/LR_DAY25-QCed.rds"
lr.day55.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/LR_DAY55-QCed.rds"

integrated.lr.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/INTEGRATED_LR_days25+55.rds"
integrated.sr.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/INTEGRATED_SR_days25_55.rds"
integrated.lr.REGRESSED.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/INTEGRATED_LR_days25+55-REGRESSED.rds"
integrated.sr.REGRESSED.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/INTEGRATED_SR_days25_55-REGRESSED.rds"

### DEFINE FUNCTIONS----------------------------------------------------------------------
auto_detect_tissue_type <- function(path_to_db_file, seuratObject, scaled, assay = "RNA", ...){
  ## THIS FUNCTION IS ADAPTED FROM THE SCTYPE GITHUB : "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R"
  ## A try-catch block was added to ensure the function continues in case an error is thrown in the for loop
  
  # get all tissue types in DB
  db_read <- openxlsx::read.xlsx(path_to_db_file)
  tissues_ <- unique(db_read$tissueType)
  result_ <- c()
  
  for(tissue in tissues_){
    print(paste0("Checking...", tissue))
    
    # Prepare gene sets
    gs_list <- gene_sets_prepare(path_to_db_file, tissue)
    
    # Prepare obj
    if(scaled){
      obj <- as.matrix(seuratObject[[assay]]@scale.data)
    } else {
      obj <- as.matrix(seuratObject[[assay]]@counts)
    }
    
    # Use tryCatch to handle errors and warnings
    tryCatch({
      es.max <- sctype_score(scRNAseqData = obj, scaled = scaled, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                             marker_sensitivity = gs_list$marker_sensitivity, verbose = FALSE)
      
      cL_resutls <- do.call("rbind", lapply(unique(seuratObject@meta.data$seurat_clusters), function(cl){
        es.max.cl <- sort(rowSums(es.max[, rownames(seuratObject@meta.data[seuratObject@meta.data$seurat_clusters==cl, ])]), decreasing = TRUE)
        head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
      }))
      
      dt_out <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1)
      
      # Return mean score for tissue
      result_ <- rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
    }, error = function(e) {
      # Handle error, e.g., print a message or set result_ to NULL
      print(paste0("An error occurred for tissue: ", tissue))
      print(e)
    }, warning = function(w) {
      # Handle warning, e.g., print a message
      print(paste0("A warning occurred for tissue: ", tissue))
      print(w)
    })
  }
  
  # Order by mean score
  result_ <- result_[order(-result_$score),]
  
  # Plot
  barplot(height = result_$score, names = result_$tissue, col = rgb(0.8, 0.1, 0.1, 0.6),
          xlab = "Tissue", ylab = "Summary score", main = "The higher summary score, the more likely tissue type is")
  
  result_
}

find_cellTypes <- function(seuratobj.path, fig_name = '', usePanglaoDB = FALSE, longRead = FALSE){
  # Load the Seurat object (from script 2) using its RDS file----------------
  seurat.obj <- readRDS(file = seuratobj.path)
  
  # We need to convert gene ids to their names for this to work, this needs to be done manually if working
  # with FLAMES
  if(longRead == TRUE) {
    # Assuming rownames contains the gene IDs
    gene_ids <- rownames(seurat.obj[["RNA"]]@scale.data)
    gene_names <- ListofGeneIDstoNames(gene_ids, gtfFilePath)
    # Set the gene names as the new row names in the RNA assay's scale.data matrix
    rownames(seurat.obj[["RNA"]]@scale.data) <- gene_names
  }
  
  # Load the desired marker database file based on usePanglaoDB parameter
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue = "Brain"
  if(usePanglaoDB == TRUE) {
    db_ = "/data/gpfs/projects/punim0646/manveer/panglao-sctypeFormat.xlsx"
  }
#  View(seurat.obj[["RNA"]]@scale.data)
  # Run the function that guesses the tissue type, we can only do this with the scType database, as it contains markers from other tissues
  if(!usePanglaoDB) {
    tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = seurat.obj, scaled = TRUE, assay = "RNA")
    tissue_guess
    tissue.guess.plt <- ggplot(tissue_guess, aes(x = tissue_guess$tissue, y = tissue_guess$score, fill = tissue)) +
      geom_bar(stat = "identity", alpha = 0.7) +
      labs(x = "Tissue", y = "Summary score (scType Database)", title = "The higher summary score, the more likely tissue type is") +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Set size and rotation for x-axis labels
        axis.text.y = element_text(size = 12),  # Set size for y-axis labels
        axis.title = element_text(size = 14),   # Set size for axis titles
        plot.title = element_text(size = 16)    # Set size for the plot title
      )
    tissue.guess.plt
  }
  
  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  gs_list
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seurat.obj[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # merge by cluster-----------------------------------------------------
  cL_results = do.call("rbind", lapply(unique(seurat.obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat.obj@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  #View(sctype_scores)
  #View(cL_results)
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  # View UMAP with cell type labels
  seurat.obj@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat.obj@meta.data$customclassif[seurat.obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  labeledUMAP <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') +
    labs(title = fig_name) +
    ggmin::theme_powerpoint() +
    theme(plot.title = element_text(size = 20))
  
  baseUMAP <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, repel = TRUE) + 
    labs(title = "Unsupervised Clustering", color = "cluster \n(from PCA)") +
    ggmin::theme_powerpoint() +
    theme(plot.title = element_text(size = 20))
  
  labeledUMAP | baseUMAP
  
  ### CREATE A BUBBLE PLOT SHOWING THE OTHER CELL TYPES THAT SCTYPE CONSIDERED--------
  # Code used here is adapted from : https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md
  
  # Prepare the edges
  cL_results = cL_results[order(cL_results$cluster),]
  edges = cL_results
  #View(cL_results)
  edges$type = paste0(edges$type, "_", edges$cluster)
  edges$cluster = paste0("cluster ", edges$cluster)
  edges = edges[, c("cluster", "type")]
  colnames(edges) = c("from", "to")
  rownames(edges) = NULL
  
  # Prepare nodes
  nodes_lvl1 = sctype_scores[, c("cluster", "ncells")]
  nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster)
  nodes_lvl1$Colour = "#f1f1ef"
  nodes_lvl1$ord = 1
  nodes_lvl1$realname = nodes_lvl1$cluster
  nodes_lvl1 = as.data.frame(nodes_lvl1)
  nodes_lvl2 = data.frame(cluster = character(),
                          ncells = numeric(),
                          Colour = character(),
                          ord = numeric(),
                          realname = character(),
                          stringsAsFactors = FALSE)
  
  ccolss = c("#5f75ae", "#92bbb8", "#64a841", "#e5486e", "#de8e06", "#eccf5a", "#b5aa0f", "#e4b680", "#7ba39d", "#b15928", "#ffff99", "#6a3d9a", "#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
  
  for (i in 1:length(unique(cL_results$cluster))) {
    dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]
    nodes_tmp = data.frame(cluster = paste0(dt_tmp$type, "_", dt_tmp$cluster),
                           ncells = dt_tmp$scores,
                           Colour = ccolss[i],
                           ord = 2,
                           realname = dt_tmp$type,
                           stringsAsFactors = FALSE)
    nodes_lvl2 = rbind(nodes_lvl2, nodes_tmp)
  }
  
  nodes = rbind(nodes_lvl1, nodes_lvl2)
  nodes$ncells[nodes$ncells < 1] = 1
  files_db = openxlsx::read.xlsx(db_)[, c("cellName", "shortName")]
  files_db = unique(files_db)
  nodes = merge(nodes, files_db, all.x = TRUE, all.y = FALSE, by.x = "realname", by.y = "cellName", sort = FALSE)
  nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]
  nodes = nodes[, c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  
  # Identify duplicate rows based on the "cluster" column
  duplicated_rows <- duplicated(nodes$cluster)
  # Filter the dataframe to keep only the first occurrence of each unique value
  nodes <- nodes[!duplicated_rows, ]
  mygraph <- graph_from_data_frame(edges, vertices = nodes)
  
  
  # Make the bubble graph
  gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
    geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
    geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
    theme_void() + 
    geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#000000"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) +
    geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
  
  bubble.plt <- scater::multiplot(DimPlot(seurat.obj, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss) + 
                                    labs(title = fig_name) + ggmin::theme_powerpoint(), gggr, cols = 2)
  
  
  ### CREATE PIE CHARTS SHOWING THE SAME INFORMATION AS THE BUBBLE PLOT----------------------
  # Create a new dataframe by filtering out rows with ord equal to 2, then remove unwanted cols
  filtered_nodes <- subset(nodes, ord != 1)
  filtered_nodes <- subset(filtered_nodes, select = -c(4, 5))
  
  # Extract the number after "_" and replace the string with the extracted number in the cluster column
  filtered_nodes$cluster <- str_replace(filtered_nodes$cluster, ".*_([0-9]+)", "Cluster \\1")
  
  # Group the dataframe by 'cluster' column
  grouped_df <- aggregate(ncells ~ realname + cluster, data = filtered_nodes, FUN = sum)
  
  # Calculate proportions within each cluster group
  grouped_df <- filtered_nodes %>%
    group_by(cluster) %>%
    mutate(prop = ncells / sum(ncells)) %>%
    ungroup() %>% 
    arrange(cluster)
  
  # Creat a list to store each pie chart generated
  pieChartList <- list()
  # Create a separate pie chart for each unique cluster group
  for (cluster in unique(grouped_df$cluster)) {
    # Subset the data for the current cluster
    cluster_data <- grouped_df[grouped_df$cluster == cluster, ]
    
    # Create the pie chart
    pie_chart <- ggplot(cluster_data, aes(x = "", y = prop, fill = realname)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      labs(title = paste(cluster)) +
      theme_void() +
      theme(legend.position = "right") +
      labs(fill = "Labels Considered")  # Rename the legend
    
    # Add the percentages to the pie chart with white color, rounded to 1 decimal place
    pie_chart <- pie_chart + 
      geom_text(aes(label = paste0(round(prop * 100, 1), "%"), x = 1.5), 
                position = position_stack(vjust = 0.5), size = 4)
    
    # Append the pie chart to the list
    pieChartList <- c(pieChartList, list(pie_chart))
  }
  pieChartList
  
  ## If working with long read data : We need replace row.names of scale.data slot back to their original gene ids
  if(longRead == TRUE) {
    # Set the row names in the RNA assay's scale.data slot back to gene ids
    rownames(seurat.obj[["RNA"]]@scale.data) <- gene_ids
  }
  
  ### Organize what the function returns------------------
  if(usePanglaoDB == TRUE){
    outputs <- list(baseUMAP, labeledUMAP, bubble.plt, pieChartList, seurat.obj)
  }else{
    outputs <- list(baseUMAP, labeledUMAP, bubble.plt, pieChartList, seurat.obj, tissue.guess.plt)
  }
  return(outputs)
}

generate_PDFs <- function(cellTypes.output.sctype, cellTypes.output.panglao, filePrefix = '', sampleDescription = ''){
  # Create a PDF to see bubble plots
  fileName_bubble <- paste0("2-",filePrefix,"-bubbles+tissueGuess-Script3.pdf")
  pdf(fileName_bubble, width = 20, height = 10)
  # Print the tissue guess bar chart
  grid.draw(cellTypes.output.sctype[[6]])
  grid.newpage()
  # Print the bubble plot for labels generated using scType database
  grid.draw(cellTypes.output.sctype[[3]])
  grid.newpage()
  # Print the bubble plot for labels generated using panglao database
  grid.draw(cellTypes.output.panglao[[3]])
  dev.off()
  
  # Open a second PDF file for exporting pie charts and UMAPs
  fileName <- paste0("1-",filePrefix, "-Script3.pdf")
  pdf(fileName, width = 20, height = 15)  # Adjust width and height as needed
  
  # Add UMAPS for comparison
  UMAP.layout <- grid.arrange(cellTypes.output.sctype[[2]], cellTypes.output.panglao[[2]],
                              nrow = 2, ncol = 2,
                              top = textGrob(paste0("SCRIPT 3 : LABELLING CELL TYPES\n",sampleDescription,"\n\n"),
                                             gp = gpar(fontsize = 25)))
  grid.draw(UMAP.layout)
  # Add Bubble Plots for scType DB labels + UMAP side to side comparison
  #grid.draw(bubble.plt1)
  UMAP.layout1 <- grid.arrange(cellTypes.output.sctype[[1]], cellTypes.output.sctype[[2]],
                               nrow = 2, ncol = 2,
                               top = textGrob("scType Database Labels versus Unsupervised Clustering\n\n",
                                              gp = gpar(fontsize = 25)))
  grid.draw(UMAP.layout1)
  
  # Add pie charts for scTypeDB labels
  piechart.layout1 <- grid.arrange(grobs = cellTypes.output.sctype[[4]], nrow = 5, ncol = 3,
                                   top = textGrob(paste0(sampleDescription,"\n(scType DataBase)\nPie Charts denoting cell label confidence scores for each cluster\n"),
                                                  gp = gpar(fontsize = 20)))
  grid.draw(piechart.layout1)
  # Add Bubble Plots for Panglao DB labels + UMAP side to side comparison
  #grid.draw(bubble.plt2)
  UMAP.layout2 <- grid.arrange(cellTypes.output.panglao[[1]], cellTypes.output.panglao[[2]],
                               nrow = 2, ncol = 2,
                               top = textGrob("Panglao Database Labels versus Unsupervised Clustering\n\n",
                                              gp = gpar(fontsize = 25)))
  grid.draw(UMAP.layout2)
  
  ## Add pie charts for panglaoDB labels to the PDF
  piechart.layout2 <- grid.arrange(grobs = cellTypes.output.panglao[[4]], nrow = 5, ncol = 3,
                                   top = textGrob(paste0(sampleDescription, "\n(Panglao DataBase)\nPie Charts denoting cell label confidence scores for each cluster\n"),
                                                  gp = gpar(fontsize = 20)))
  grid.draw(piechart.layout2)
  # Close the PDF device
  dev.off()
}

#### Run the find_cellType function on each seurat object with both scType and Panglao databases, 
#### then generate pdfs to display results------------------------------------------------------------------------------
integrated.sr.sctypeDB.REGRESSED <- find_cellTypes(integrated.sr.REGRESSED.filepath, "Integrated Short Reads (REGRESSED) - scTypeDB")
integrated.sr.panglaoDB.REGRESSED <- find_cellTypes(integrated.sr.REGRESSED.filepath, "Integrated Short Reads (REGRESSED) - PanglaoDB", usePanglaoDB = TRUE)
generate_PDFs(integrated.sr.sctypeDB.REGRESSED, integrated.sr.panglaoDB.REGRESSED, 
              filePrefix = "Integrated_SR_REGRESSED", sampleDescription = "INTEGRATED SHORT READS (REGRESSED)")

integrated.lr.sctypeDB.REGRESSED <- find_cellTypes(integrated.lr.REGRESSED.filepath, "Integrated Long Reads (REGRESSED) - scTypeDB", longRead = TRUE)
integrated.lr.panglaoDB.REGRESSED <- find_cellTypes(integrated.lr.REGRESSED.filepath, "Integrated Long Reads (REGRESSED) - PanglaoDB", usePanglaoDB = TRUE,
                                          longRead = TRUE)
generate_PDFs(integrated.lr.sctypeDB.REGRESSED, integrated.lr.panglaoDB.REGRESSED, filePrefix = "Integrated_LR_REGRESSED", sampleDescription = "INTEGRATED LONG READS")

sr.day25.sctypeDB <- find_cellTypes(sr.day25.filepath, "Day 25 (Short Read) - scTypeDB")
sr.day25.panglaoDB <- find_cellTypes(sr.day25.filepath, "Day 25 (Short Read) - PanglaoDB", usePanglaoDB = TRUE)
#generate_PDFs(sr.day25.sctypeDB, sr.day25.panglaoDB, filePrefix = "SR-day25", sampleDescription = "DAY 25 (SHORT READ)")

sr.day55.sctypeDB <- find_cellTypes(sr.day55.filepath, "Day 55 (Short Read) - scTypeDB")
sr.day55.panglaoDB <- find_cellTypes(sr.day55.filepath, "Day 55 (Short Read) - PanglaoDB", usePanglaoDB = TRUE)
#generate_PDFs(sr.day55.sctypeDB, sr.day55.panglaoDB, filePrefix = "SR-day55", sampleDescription = "DAY 55 (SHORT READ)")

q20.day25.sctypeDB <- find_cellTypes(q20.day25.filepath, "Day 25 (Q20 Long Read) - scTypeDB", longRead = TRUE)
q20.day25.panglaoDB <- find_cellTypes(q20.day25.filepath, "Day 25 (Q20 Long Read) - PanglaoDB", usePanglaoDB = TRUE, longRead = TRUE)
#generate_PDFs(q20.day25.sctypeDB, q20.day25.panglaoDB, filePrefix = "Q20-day25", sampleDescription = "DAY 25 (LONG READ Q20)")

lr.day25.sctypeDB <- find_cellTypes(lr.day25.filepath, "Day 25 (Long Read) - scTypeDB", longRead = TRUE)
lr.day25.panglaoDB <- find_cellTypes(lr.day25.filepath, "Day 25 (Long Read) - PanglaoBD", usePanglaoDB = TRUE, longRead = TRUE)
#generate_PDFs(lr.day25.sctypeDB, lr.day25.panglaoDB, filePrefix = "LR-DAY25", sampleDescription = "Day 25 (LONG READ PROMETHION)")

lr.day55.sctypeDB <- find_cellTypes(lr.day55.filepath, "Day 55 (Long Read) - scTypeDB", longRead = TRUE)
lr.day55.panglaoDB <- find_cellTypes(lr.day55.filepath, "Day 55 (Long Read) - PanglaoBD", usePanglaoDB = TRUE, longRead = TRUE)
#generate_PDFs(lr.day55.sctypeDB, lr.day55.panglaoDB, filePrefix = "LR-DAY55", sampleDescription = "Day 55 (LONG READ PROMETHION)")


### Extract seurat objects with their cell labels - HEATMAP-------------------------------------------------------------------------
seurat.obj.lr.day55.sctypeDB <- lr.day55.sctypeDB[[5]]
View(seurat.obj.lr.day55.sctypeDB@meta.data)

# find markers for every cluster compared to all remaining cells, report only the positive ones
lr.day55.markers <- FindAllMarkers(seurat.obj.lr.day55.sctypeDB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Group markers by cluster, then keep the top 3 features with the highest log fold change for each cluster
top3 <- lr.day55.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

DoHeatmap(seurat.obj.lr.day55.sctypeDB, features = top3$gene)

#Idents(seurat.obj.lr.day55.sctypeDB) <- "customclassif"
#lr.day55.markers <- FindAllMarkers(seurat.obj.lr.day55.sctypeDB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top3 <- lr.day55.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 3, order_by = avg_log2FC)

#DoHeatmap(seurat.obj.lr.day55.sctypeDB, features = top3$gene)









### Investigate canonical marker genes of expected cell types------------------------------

## Extract seurat objects with their labels :----------------------------------------------
seurat.obj.sr.day55.sctypeDB <- sr.day55.sctypeDB[[5]]
seurat.obj.sr.day25.sctypeDB <- sr.day25.sctypeDB[[5]]
seurat.obj.lr.day55.sctypeDB <- lr.day55.sctypeDB[[5]]
seurat.obj.lr.day25.sctypeDB <- lr.day25.sctypeDB[[5]]
seurat.obj.q20.day25.sctypeDB <- q20.day25.sctypeDB[[5]]

seurat.obj.integratedSR <- integrated.sr.sctypeDB[[5]]
seurat.obj.integratedLR <- integrated.lr.sctypeDB[[5]]
seurat.obj.integratedSR.REGRESSED <- integrated.sr.sctypeDB.REGRESSED[[5]]
seurat.obj.integratedLR.REGRESSED <- integrated.lr.sctypeDB.REGRESSED[[5]]

sr.day55.sctypeDB[[1]] | sr.day55.sctypeDB[[2]]

### Define functions for investigating the expression of marker gene sets, and generating a pdf-------------------
FeatureSetPlot <- function(seurat.obj, gene.marker.list, plot.name){
  # https://github.com/satijalab/seurat/issues/3521
  # https://github.com/satijalab/seurat/issues/4544
  # https://satijalab.org/seurat/reference/addmodulescore
  seurat.obj <- AddModuleScore(object = seurat.obj,
                               features = gene.marker.list,
                               name = plot.name)
  plt <- FeaturePlot(object = seurat.obj, features = paste0(plot.name, "1"))
  
  plt
}
generateGeneMarkerPlots <- function(seurat.obj, sampleDesc = "", longRead = FALSE){
  
  Excitatory_deep_layer <- list(c("SOX5", "NR4A2", "CRYM", "TBR1", "FOXP2"))
  Deep_layer_markers <- list(c("RORB", "FOXP2", "CRYM", "TBR1", "NR4A2"))
  Maturing_excitatory_markers <- list(c("SATB2", "STMN2", "NEUROD6"))
  #Inhibitory_interneuron_markers <- list(c("DLX1", "DLX2", "LHX6", "DLX5")) # a single cell transcriptomic atlas of human neocortical development paper
  
  Canonical_Radial_Glial_markers <- list(c("VIM", "HES1", "PAX6", "SOX2"))
  GABAergic.markers <- list(c("GAD1", "GAD2", "SST", "DLX1", "DLX2", "DLX5"))
  Astrocyte_markers <- list(c("S100B", "GJA1", "AQP4", "GFAP"))
  Glutamatergic.markers <- list(c("GRIN2A", "GRIN2B", "GRIN1", "GRIA1", "GRIA2",
                                  "GRIA3", "GRIA4", "SLC17A7", "NEUROD6", "SLC1A2"))
 
  # We need to convert our gene markers from name format to their respective symbols for long-reads
  if(longRead == TRUE) {
    Canonical_Radial_Glial_markers <- list(ListofGeneNamestoIDs(Canonical_Radial_Glial_markers[[1]], gtfFilePath))
    GABAergic.markers <- list(ListofGeneNamestoIDs(GABAergic.markers[[1]], gtfFilePath))
    Astrocyte_markers <- list(ListofGeneNamestoIDs(Astrocyte_markers[[1]], gtfFilePath))
    Glutamatergic.markers <- list(ListofGeneNamestoIDs(Glutamatergic.markers[[1]], gtfFilePath))
  }
  
  Radial.Glia.plt <- FeatureSetPlot(seurat.obj, Canonical_Radial_Glial_markers, "Radial_glia_genes") +
    labs(title = "Radial Glia Gene Markers: \nVIM, HES1, PAX6, SOX2")
  
  GABAergic.plt <- FeatureSetPlot(seurat.obj, GABAergic.markers, "GABAergic_genes") +
    labs(title = "GABAergic Gene Markers: \nGAD1, GAD2, SST,\nDLX1, DLX2, DLX5")
  
  Astrocyte.plt <- FeatureSetPlot(seurat.obj, Astrocyte_markers, "Astrocyte_genes") +
    labs(title = "Astrocyte Gene Markers:\nS100B, GJA1, AQP4, GFAP")
  
  Glutamatergic.plt <- FeatureSetPlot(seurat.obj, Glutamatergic.markers, "Glutamatergic_genes") +
    labs(title = "Glutamatergic Gene Markers:\nGRIN1, GRIN2A, GRIN2B,\nGRIA1, GRIA2, GRIA3, GRIA4,\nNEUROD6, SLC1A7, SLC1A2")
  
  if(longRead == TRUE) {
    EOMES.plt <- FeaturePlot(seurat.obj, features = "ENSG00000163508.13") +
      labs(title = "EOMES\n(Intermediate Progenitor Cell Marker)")
  } else {
    EOMES.plt <- FeaturePlot(seurat.obj, features = "EOMES") +
      labs(title = "EOMES\n(Intermediate Progenitor Cell Marker)")
  }
  
  marker.plots.layout <- grid.arrange(Radial.Glia.plt, Astrocyte.plt,
                                      Glutamatergic.plt, GABAergic.plt,
                                      EOMES.plt,
                                      nrow = 3, ncol = 2,
                                      top = textGrob(paste0("Cell Type Gene Markers : ", sampleDesc),
                                                     gp = gpar(fontsize = 20)))
  
  pdf.figs.title <- paste0("3-",sampleDesc, "-Gene-Markers-script3.pdf")
  ggsave(pdf.figs.title, marker.plots.layout, width = 15, height = 20)
}

## Missing genes for long read day 55-----------------------------------
missing_genes <- c("ENSG00000128683.14",
                   "ENSG00000171885.18",
                   "ENSG00000183454.18", 
                   "ENSG00000176884.17")
missing_genes <- ListofGeneIDstoNames(missing_genes, gtfFilePath)
missing_genes # LR didn't detect genes : GAD1, AQP4, GRIN2A, GRIN1 for day 55

## Generate gene marker plots for DAY 55 LONG READ AND SHORT READS (NOT WORKING FOR DAY 25 - MISSING GENE PROBLEM)
generateGeneMarkerPlots(seurat.obj.integratedSR.REGRESSED, "Integrated-SR-REGRESSED")
generateGeneMarkerPlots(seurat.obj.integratedLR.REGRESSED, "Integrated-LR-REGRESSED", longRead = TRUE)

generateGeneMarkerPlots(seurat.obj.lr.day55.sctypeDB, "LR-day55", longRead = TRUE)
generateGeneMarkerPlots(seurat.obj.sr.day55.sctypeDB, "SR-day55")
#generateGeneMarkerPlots(seurat.obj.sr.day25.sctypeDB, "SR-day25")
#generateGeneMarkerPlots(seurat.obj.lr.day25.sctypeDB, "LR-day25", longRead = TRUE)
#generateGeneMarkerPlots(seurat.obj.q20.day25.sctypeDB, "q20-day25", longRead = TRUE)
lr.day55.sctypeDB[[2]] | lr.day55.panglaoDB[[2]]

FeaturePlot(seurat.obj.sr.day55.sctypeDB, features = "NEUROG2")
FeaturePlot(seurat.obj.sr.day55.sctypeDB, features = "PAX6")
sr.day25.sctypeDB[[2]]
FeaturePlot(seurat.obj.sr.day55.sctypeDB, features = "PDGFA")

# EOMES - marker of intermediate progenitor cells
FeaturePlot(seurat.obj.integratedLR, features = "ENSG00000163508.13") +
  labs(title = "EOMES - Integrated Long Reads")
FeaturePlot(seurat.obj.integratedSR, features = "EOMES") +
  labs(title = "EOMES - Integrated Short Reads")
ListofGeneNamestoIDs(list("EOMES"), gtfFilePath)

# VERSION 2 DOESNT WORK BECAUSE NOT ALL GENES ARE FOUND WITHIN LONG READ SAMPLES---------------------------
# Version 2 that uses PercentangeFeatureSet() instead of AddModuleScore()
generateGeneMarkerPlots.v2 <- function(seurat.obj, sampleDesc = "", longRead = FALSE){
  seurat.obj <- seurat.obj.lr.day55.sctypeDB
  longRead <- TRUE
  
  Excitatory_deep_layer <- c("SOX5", "NR4A2", "CRYM", "TBR1", "FOXP2")
  Deep_layer_markers <- c("RORB", "FOXP2", "CRYM", "TBR1", "NR4A2")
  Maturing_excitatory_markers <- c("SATB2", "STMN2", "NEUROD6")
  #Inhibitory_interneuron_markers <- c("DLX1", "DLX2", "LHX6", "DLX5") # a single cell transcriptomic atlas of human neocortical development paper
  
  Canonical.Radial.Glial.markers <- c("VIM", "HES1", "PAX6", "SOX2")
  GABAergic.markers <- c("GAD1", "GAD2", "SST", "DLX1", "DLX2", "DLX5")
  Astrocyte.markers <- c("S100B", "GJA1", "AQP4", "GFAP")
  Glutamatergic.markers <- c("GRIN2A", "GRIN2B", "GRIN1", "GRIA1", "GRIA2",
                                  "GRIA3", "GRIA4", "SLC17A7", "NEUROD6", "SLC1A2")
  
  # We need to convert our gene markers from name format to their respective symbols for long-reads
  if(longRead == TRUE) {
    Canonical.Radial.Glial.markers <- ListofGeneNamestoIDs(Canonical.Radial.Glial.markers, gtfFilePath)
    GABAergic.markers <- ListofGeneNamestoIDs(GABAergic.markers, gtfFilePath)
    Astrocyte.markers <- ListofGeneNamestoIDs(Astrocyte.markers, gtfFilePath)
    Glutamatergic.markers <- ListofGeneNamestoIDs(Glutamatergic.markers, gtfFilePath)
  }
  Canonical.Radial.Glial.markers
  FeaturePlot(seurat.obj, feature = GABAergic.markers[7])
  
  ## Radial Glia markers percentage feature set plot
  seurat.obj[["Radial.Glia.Markers"]] <- PercentageFeatureSet(seurat.obj,
                                                              features = Canonical.Radial.Glial.markers)
  Radial.Glia.plt <- FeaturePlot(seurat.obj, reduction = "umap", features = 'Radial.Glia.Markers')  +
    labs(title = "Radial Glia Markers: \nVIM, HES1, PAX6, SOX2")
  
  ## GABAergic markers percentage feature set plot
  seurat.obj[["GABAergic.Neuron.Markers"]] <- PercentageFeatureSet(seurat.obj,
                                                              features = GABAergic.markers[2:3])
  GABAergic.plt <- FeaturePlot(seurat.obj, reduction = "umap", features = 'GABAergic.Neuron.Markers')  +
    labs(title = "GABAergic Markers: \nGAD1, GAD2, SST,\nDLX1, DLX2, DLX5")
  
  ## Astrocyte markers percentage feature set plot
  seurat.obj[["Astrocyte Markers"]] <- PercentageFeatureSet(seurat.obj,
                                                                   features = Astrocyte.markers)
  Astrocyte.plt <- FeaturePlot(seurat.obj, reduction = "umap", features = 'Astrocyte Markers')  +
    labs(title = "Astrocyte Markers:\nS100B, GJA1, AQP4, GFAP")
  
  ## Glutamatergic markers percentage feature set plot
  seurat.obj[["Glutamatergic Neuron Markers"]] <- PercentageFeatureSet(seurat.obj,
                                                            features = Glutamatergic.markers)
  Glutamatergic.plt <- FeaturePlot(seurat.obj, reduction = "umap", features = 'Glutamatergic Neuron Markers')  +
    labs(title = "Glutamatergic Markers:\nGRIN1, GRIN2A, GRIN2B,\nGRIA1, GRIA2, GRIA3, GRIA4,\nNEUROD6, SLC1A7, SLC1A2P")

  marker.plots.layout <- grid.arrange(Radial.Glia.plt, Astrocyte.plt,
                                      Glutamatergic.plt, GABAergic.plt,
                                      nrow = 2, ncol = 2,
                                      top = textGrob(paste0("Cell Type Gene Markers : ", sampleDesc),
                                                     gp = gpar(fontsize = 20)))
  
  pdf.figs.title <- paste0("3-",sampleDesc, "-Gene-Markers-script3-v2.pdf")
  ggsave(pdf.figs.title, marker.plots.layout, width = 15, height = 15)
}
generateGeneMarkerPlots.v2(seurat.obj.lr.day55.sctypeDB, "LR-day55")

Neuronal.differentiation.markers <- list(c("ASCL1", "CD24", "CDH2", "EPCAM",
                                           "ITGB1", "L1CAM", "MAP2",
                                           "MCAM", "NEUROG2", "SATB2",
                                           "SOX5", "TUBB3"))
neuron.diff.marker.plt <- FeatureSetPlot(seurat.obj.sr.day55.sctypeDB, Neuronal.differentiation.markers, "neuron_diff") +
  labs(title = "Neuronal Differentiation Markers")
neuron.diff.marker.plt

neuronal.progenitor.cell.markers <- list(c("CUX1", "CUX2", "EMX2", "EOMES", "FEZF2",
                                          "LMO4", "SVET1", "TBR1", "BCL11B"))
neural.progenitor.marker.plt <- FeatureSetPlot(seurat.obj.sr.day55.sctypeDB, neuronal.progenitor.cell.markers, "neuron.progenitor") +
  labs(title = "Neuronal Progenitor Cell Differentiation Markers")
neural.progenitor.marker.plt




#### TINKER WITH LABELS - BASED ON INFORMATION THUS FAR
# Change long read cell labels
seurat.obj.integratedLR.REGRESSED$customclassif[seurat.obj.integratedLR.REGRESSED$seurat_clusters == "6"] <- "Neural stem cells"
seurat.obj.integratedLR.REGRESSED$customclassif[seurat.obj.integratedLR.REGRESSED$seurat_clusters == "9"] <- "Neural stem cells"
seurat.obj.integratedLR.REGRESSED$customclassif[seurat.obj.integratedLR.REGRESSED$seurat_clusters == "7"] <- "Glutamatergic neurons"
seurat.obj.integratedLR.REGRESSED$customclassif[seurat.obj.integratedLR.REGRESSED$seurat_clusters == "0"] <- "Glutamatergic neurons"
seurat.obj.integratedLR.REGRESSED$customclassif[seurat.obj.integratedLR.REGRESSED$seurat_clusters == "4"] <- "Radial glial cells"
seurat.obj.integratedLR.REGRESSED$customclassif[seurat.obj.integratedLR.REGRESSED$seurat_clusters == "3"] <- "Radial glial cells"

# Change short read cell labels
seurat.obj.integratedSR.REGRESSED$customclassif[seurat.obj.integratedSR.REGRESSED$seurat_clusters == "6"] <- "Radial glial cells"
seurat.obj.integratedSR.REGRESSED$customclassif[seurat.obj.integratedSR.REGRESSED$seurat_clusters == "3"] <- "Glutamatergic neurons"

# Plot our relabeled UMAPS
LR.relabeled.UMAP <- DimPlot(seurat.obj.integratedLR.REGRESSED, reduction = 'umap', 
                         group.by = 'customclassif', label = TRUE,
                         label.size = 5, repel = TRUE) +
  labs(title = "Integrated Long Reads\n(Re-evaluated Labels)") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
LR.relabeled.UMAP | integrated.lr.sctypeDB.REGRESSED[[2]]

SR.relabeled.UMAP <- DimPlot(seurat.obj.integratedSR.REGRESSED, reduction = 'umap', 
                             group.by = 'customclassif', label = TRUE,
                             label.size = 5, repel = TRUE) +
  labs(title = "Integrated Short Reads\n(Re-evaluated Labels)") +
  ggmin::theme_powerpoint() +
  theme(title = element_text(size = 13))
SR.relabeled.UMAP | integrated.sr.sctypeDB.REGRESSED[[2]]


saveRDS(seurat.obj.integratedLR.REGRESSED, file = "Relabeled_LR_days25+55-script3.rds")
saveRDS(seurat.obj.integratedSR.REGRESSED, file = "Relabeled_SR_days25+55-script3.rds")
