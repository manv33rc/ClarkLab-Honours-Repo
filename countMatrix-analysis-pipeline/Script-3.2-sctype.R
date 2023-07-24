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

set.seed(4242)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

sr.day25.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/filtered_sr_Day25-QCed.rds"
sr.day55.filepath <- "/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/sr_Day55_FIXED-QCed.rds"

find_cellTypes <- function(seuratobj.path, fig_name = '', usePangoDB = FALSE){
  # Load the Seurat object (from script 2) using its RDS file----------------
  seurat.obj <- readRDS(file = seuratobj.path)
  
  # DB file
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue = "Brain"
  if(usePangoDB == TRUE) {
    db_ = "/data/gpfs/projects/punim0646/manveer/panglao-sctypeFormat.xlsx"
  }
  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  
  
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
  #View(nodes)
  
  mygraph <- graph_from_data_frame(edges, vertices = nodes)
  
  
  # Make the bubble graph
  gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
    geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
    geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
    theme_void() + 
    geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) + 
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

  
  ## Display all the pie charts generated in a grid layout
#  n_pie_charts <- length(pieChartList)
#  n_cols <- min(3, n_pie_charts) + 1
#  n_rows <- ceiling(n_pie_charts / n_cols)
  
#  pieChart.layout <- grid.arrange(grobs = pieChartList, ncol = n_cols, 
#                                  layout_matrix = matrix(seq_len(n_pie_charts), 
#                                                         nrow = n_rows, ncol = n_cols, 
#                                                         byrow = TRUE),
#                                  top = textGrob(paste(fig_name, " Pie Charts"),
#                                                 gp = gpar(fontsize = 20)))
  
  ### Organize what the function returns------------------
  outputs <- list(baseUMAP, labeledUMAP, bubble.plt, pieChartList, seurat.obj)
  return(outputs)
}

#### Run the find_cellType functions on each seurat object-----------------------
day25.sctypeDB <- find_cellTypes(sr.day25.filepath, "Day 25 (Short Read) - scTypeDB")
day25.panglaoDB <- find_cellTypes(sr.day25.filepath, "Day 25 (Short Read) - PanglaoDB", usePangoDB = TRUE)

day55.sctypeDB <- find_cellTypes(sr.day55.filepath, "Day 55 (Short Read) - scTypeDB")
day55.panglaoDB <- find_cellTypes(sr.day55.filepath, "Day 55 (Short Read) - PanglaoDB", usePangoDB = TRUE)

### View UMAPS (labeled and base) from panglao and sctype's marker databases for SR day 25 and day 55------------
day25.sctypeDB[[2]] | day25.panglaoDB[[2]] | day25.sctypeDB[[1]]
day25.sctypeDB[[2]] | day25.panglaoDB[[2]]

### View bubble charts from both panglaoDB and sctype for SR day 25----------------------------
bubble.plt1 <- day25.sctypeDB[[3]]
grid.newpage()
grid.draw(bubble.plt1)
bubble.plt2 <- day25.panglaoDB[[3]]
grid.newpage()
grid.draw(bubble.plt2)

# Create a grid layout with two columns
bubble_grid <- grid.arrange(bubble.plt1, bubble.plt2, ncol = 2)
# Display the grid layout
grid.newpage()
grid.draw(bubble_grid)

### View bubble charts from both panglaoDB and sctype for SR day 55----------------------------
bubble.plt1 <- day55.sctypeDB[[3]]
grid.newpage()
grid.draw(bubble.plt1)
bubble.plt2 <- day55.panglaoDB[[3]]
grid.newpage()
grid.draw(bubble.plt2)
# Create a grid layout with two columns
bubble_grid <- grid.arrange(bubble.plt1, bubble.plt2, ncol = 2)
# Display the grid layout
grid.newpage()
grid.draw(bubble_grid)

### View the pie charts for SR day 25----------------------------------
piecharts1 <- day25.sctypeDB[[4]]
length(piecharts1)
grid.newpage()
day25.sctypeDB[[1]] | day25.sctypeDB[[2]]
grid.arrange(grobs = piecharts1, nrow = 3, ncol = 2,
             top = textGrob("SHORT READ - Day 25 \n(scType DataBase)\nPie Charts denoting cell label confidence scores for each cluster\n",
                            gp = gpar(fontsize = 12, fontface = "bold")))

piecharts2 <- day25.panglaoDB[[4]]
grid.newpage()
day25.panglaoDB[[1]] | day25.panglaoDB[[2]]
grid.arrange(grobs = piecharts2, nrow = 3, ncol = 2,
             top = textGrob("SHORT READ - Day 25 \n(Panglao DataBase)\nPie Charts denoting cell label confidence scores for each cluster\n",
                            gp = gpar(fontsize = 12, fontface = "bold")))

### CREATE A PDF DISPLAYING OUTPUTS FOR SR DAY 55----------------------------------
# Open a PDF file for exporting
fileName <- "SR-day55-Script3.pdf"
pdf(fileName, width = 20, height = 15)  # Adjust width and height as needed

# Add UMAPS for comparison
UMAP.layout <- grid.arrange(day55.sctypeDB[[2]], day55.panglaoDB[[2]],
                            nrow = 2, ncol = 2,
                            top = textGrob("SCRIPT 3 : LABELLING CELL TYPES\nDAY 55 (SHORT READ)\n\n",
                                           gp = gpar(fontsize = 25)))
grid.draw(UMAP.layout)
# Add Bubble Plots for scType DB labels + UMAP side to side comparison
#grid.draw(bubble.plt1)
UMAP.layout1 <- grid.arrange(day55.sctypeDB[[1]], day55.sctypeDB[[2]],
                             nrow = 2, ncol = 2,
                             top = textGrob("scType Database Labels versus Unsupervised Clustering\n\n",
                                            gp = gpar(fontsize = 25)))
grid.draw(UMAP.layout1)

# Add pie charts for scTypeDB labels
piechart.layout1 <- grid.arrange(grobs = day55.sctypeDB[[4]], nrow = 4, ncol = 2,
                                 top = textGrob("SHORT READ - Day 55 \n(scType DataBase)\nPie Charts denoting cell label confidence scores for each cluster\n",
                                                gp = gpar(fontsize = 20)))
grid.draw(piechart.layout1)
# Add Bubble Plots for Panglao DB labels + UMAP side to side comparison
#grid.draw(bubble.plt2)
UMAP.layout2 <- grid.arrange(day55.panglaoDB[[1]], day55.panglaoDB[[2]],
                             nrow = 2, ncol = 2,
                             top = textGrob("Panglao Database Labels versus Unsupervised Clustering\n\n",
                                            gp = gpar(fontsize = 25)))
grid.draw(UMAP.layout2)

## Add pie charts for panglaoDB labels to the PDF
piechart.layout2 <- grid.arrange(grobs = day55.panglaoDB[[4]], nrow = 4, ncol = 2,
                                 top = textGrob("SHORT READ - Day 55 \n(Panglao DataBase)\nPie Charts denoting cell label confidence scores for each cluster\n",
                                                gp = gpar(fontsize = 20)))
grid.draw(piechart.layout2)
# Close the PDF device
dev.off()
