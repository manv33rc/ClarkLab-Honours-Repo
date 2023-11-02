# Script to find and visualize all the transcripts when provided a gene of interest

setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

set.seed(4242)

library(reticulate)
use_python("/Users/manveerchuahan/miniconda3/bin/python3.11")
library(tidyverse)
library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)

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

# Transfer cell labels from gene level object to isoform level object for long reads
integrated.lr.ISOLVL.seurat.obj@meta.data$customclassif <- integrated.lr.seurat.obj@meta.data$customclassif

DimPlot(integrated.lr.ISOLVL.seurat.obj, reduction = 'umap', group.by = 'customclassif',
        label = TRUE, label.size = 4, repel = TRUE) +
  labs(title = 'Cell Type labels for Isoform Level Long Read Data')

# Add the a column that displays the number of isoforms detected for a given gene
transcript.gene.dict <- transcript.gene.dict %>%
  mutate(isoformNumber = str_count(transcript_id, ",") + 1)
View(transcript.gene.dict)

# Define a function that returns a filtered transcript gene dictionary displaying genes with a specified isoform count
filterGeneIsoformDictionary <- function(transcript.gene.dict,
                                        target.isoformNumber = 1, exactNumber = FALSE,
                                        includesNovelIsoform = FALSE){
  if(exactNumber){
    filtered.transcript.gene.dict <- transcript.gene.dict %>% 
      filter(isoformNumber == target.isoformNumber)
  }else{
    filtered.transcript.gene.dict <- transcript.gene.dict %>%
      filter(isoformNumber >= target.isoformNumber)
  }
  
  if(includesNovelIsoform) {
    filtered.transcript.gene.dict <- filtered.transcript.gene.dict %>% 
      filter(grepl('BambuTx', transcript_id))
  }
  
  return(filtered.transcript.gene.dict)
}

# Take a look at genes with greater than 6 isoforms (and at least one novel isoform from bambu)---------------
genes.with.greater.than.6.isoforms.w.novel <- filterGeneIsoformDictionary(transcript.gene.dict = transcript.gene.dict,
                                                                          target.isoformNumber = 6,
                                                                          includesNovelIsoform = TRUE)
View(genes.with.greater.than.6.isoforms.w.novel)
# Interesting genes here : APBA2, POLB, AAAS, BIN1*, METTL23*



# Genes with known isoform switching events in the literature----------
# PTBP genes, SRRM4, GRIA2, DEF1

# Write a function that creates UMAPS and alternative visualizations for a gene of interest-------------
printIsoformsforGene <- function(lr.geneLvl.seurat.obj,
                                 lr.isoformLvl.seurat.obj,
                                 transcript.gene.dict,
                                 GOI,
                                 rowNum = 1,
                                 ridge.plt = FALSE,
                                 vln.plt = FALSE,
                                 byCellLabel = FALSE){
  # Assign base UMAP plots as a reference
  lr.celltype.geneLvl <- DimPlot(lr.geneLvl.seurat.obj, reduction = 'umap', group.by = 'customclassif',
                         label = TRUE, label.size = 4, repel = TRUE) +
    labs(title = 'Integrated Long Reads\nscType Generated Labels\n(Gene Level)')
  lr.unsupervised.geneLvl <- DimPlot(lr.geneLvl.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters',
                             label = TRUE, label.size = 5, repel = TRUE) +
    labs(title = 'Integrated Long Reads\nUnsupervised Seurat Clustering\n(Gene Level)')
  
  lr.celltype.isoLvl <- DimPlot(lr.isoformLvl.seurat.obj, reduction = 'umap', group.by = 'customclassif',
                                 label = TRUE, label.size = 4, repel = TRUE) +
    labs(title = 'Integrated Long Reads\nscType Generated Labels\n(Isoform Level)')
  lr.unsupervised.isoLvl <- DimPlot(lr.isoformLvl.seurat.obj, reduction = 'umap', group.by = 'seurat_clusters',
                                     label = TRUE, label.size = 5, repel = TRUE) +
    labs(title = 'Integrated Long Reads\nUnsupervised Seurat Clustering\n(Isoform Level)')
  
  lr.celltype.geneLvl | lr.unsupervised.geneLvl
  lr.celltype.isoLvl | lr.unsupervised.isoLvl
  
  # Create a list of the transcripts for the given gene of interest (GOI)
  GOI.row <- transcript.gene.dict[transcript.gene.dict$geneName == GOI, ]
  GOI.row <- GOI.row[rowNum, ]
  GOI.transcript.list <- strsplit(GOI.row$transcript_id, split = ',')[[1]]
  GOI.transcriptName.list <- ListofTranscriptIDstoNames(GOI.transcript.list, gtfFilePath)
  # Create an empty list to store UMAP plots of isoforms
  isoform.UMAP.plts <- list()
  iso.ridge.plts <- list()
  iso.vln.plts <- list()
  
  ## Compile UMAP plots showing the expression of all the isoforms of the GOI
  for (isoform in seq_along(GOI.transcript.list)) {
    isoformID <- GOI.transcript.list[[isoform]]
    isoformName <- GOI.transcriptName.list[[isoform]]
    
    if(is.null(isoformName)) {
      isoformName <- paste("Novel", GOI, "Isoform:", isoformID)
    }
    
    print(isoformID)
    print(isoformName)
    
    # Use tryCatch to handle errors
    tryCatch({
      # Your original code that might throw an error
      isoform.UMAP.plts[[isoform]] <- FeaturePlot(lr.isoformLvl.seurat.obj,
                                                  reduction = 'umap',
                                                  features = isoformID,
                                                  order = TRUE) +
        labs(title = paste0(isoformName, "\n(", isoformID, ")"))
      if(ridge.plt | vln.plt){
        # Set identity of seurat object to cell type labels or unsupervised cluster labels depending on what is specified
        if(byCellLabel){
          Idents(lr.isoformLvl.seurat.obj) <- lr.isoformLvl.seurat.obj$customclassif
        } else{
          Idents(lr.isoformLvl.seurat.obj) <- lr.isoformLvl.seurat.obj$seurat_clusters
        }
        
        # Generate a ridge plot or violin plot showing isoform expression distributions if specified
        if(ridge.plt){
          iso.ridge.plts[[isoform]] <- RidgePlot(lr.isoformLvl.seurat.obj,
                                                        features = isoformID,
                                                        ncol = 2) +
            labs(title = paste0(isoformName, "\n(", isoformID, ")"))
        }
        if(vln.plt){
          iso.vln.plts[[isoform]] <- VlnPlot(lr.isoformLvl.seurat.obj,
                                                        features = isoformID) +
            labs(title = paste0(isoformName, "\n(", isoformID, ")"))
        }
      }
    },
    error = function(e) {
      # What to do if an error is caught
      message("Caught an error: ", e$message)
      message("Continuing to next iteration...")
    })
  }
  
  # Create a UMAP showing the expression of the GOI (in the gene level seurat object)
  if(length(row.names(GOI.row)) > 1){
    message(paste0(length(row.names(GOI.row)), " rows found in dictionary for: ", GOI))
  }
  GOI.UMAP <- FeaturePlot(lr.geneLvl.seurat.obj,
                          reduction = 'umap',
                          features =  as.character(row.names(GOI.row)[rowNum]),
                          order = TRUE) +
    labs(title = paste0("Gene Level Expression : ", GOI, 
                       "\n(", as.character(row.names(GOI.row)[rowNum]), ")"))
  
  
  pdfFileName = paste0(GOI, "-isoformDE.pdf")
  # Open a new PDF device--------------
  pdf(file = pdfFileName, width = 11, height = 8)
  
  grid.draw(lr.celltype.geneLvl | GOI.UMAP)
  grid.draw(lr.celltype.isoLvl)
  
  for(i in seq_along(isoform.UMAP.plts)){
    print(isoform.UMAP.plts[[i]])
    if(ridge.plt){
      print(iso.ridge.plts[[i]])
    }
    if(vln.plt){
      print(iso.vln.plts[[i]])
    }
  }
  # Close the PDF device
  dev.off()
  
  message(paste0("PDF Created for ", GOI, ". ", length(isoform.UMAP.plts),  " Isoforms detected."))
  #return(isoform.UMAP.plts)
}

printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "RBFOX2")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "RBFOX3")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "RBFOX1")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "GRIA2",
                     vln.plt = TRUE,
                     byCellLabel = TRUE)
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "NIN",
                     vln.plt = TRUE,
                     ridge.plt = TRUE,
                     byCellLabel = TRUE)

printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "PTBP2")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "PTBP1")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "NOVA1")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "NUMB")

printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "METTL23")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "BIN1")
printIsoformsforGene(lr.geneLvl.seurat.obj = integrated.lr.seurat.obj,
                     lr.isoformLvl.seurat.obj = integrated.lr.ISOLVL.seurat.obj,
                     transcript.gene.dict = transcript.gene.dict,
                     GOI = "APBA2")