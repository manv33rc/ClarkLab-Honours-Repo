# Script to create a stacked bar chart to compare cell barcodes shared between long reads and short reads

setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")

set.seed(4242)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(readr)
library(RColorBrewer)


## Read in filtered seurat objects for both samples with short reads and long reads---------------
sr.day25.filepath <- "SR_DAY25-QCed.rds"
sr.day55.filepath <- "SR_DAY55-QCed.rds"
lr.day25.filepath <- "LR_DAY25-QCed.rds"
lr.day55.filepath <- "LR_DAY55-QCed.rds"

# Use filepaths to .RDS objects, so that we have the seurat objects we need to work with in this environment
short.read.day25 <- readRDS(file = sr.day25.filepath)
short.read.day55 <- readRDS(file = sr.day55.filepath)
long.read.day25 <- readRDS(file = lr.day25.filepath)
long.read.day55 <- readRDS(file = lr.day55.filepath)

# Long read barcodes
lr.D25.barcodes.path <- "/Volumes/Expansion/data1/Day-25+55_whitelists/lr-cortdiff-day25-whitelist.csv"
lr.D55.barcodes.path <- "/Volumes/Expansion/data1/Day-25+55_whitelists/lr-cortdiff-day55-whitelist.csv"


# Short read barcodes are in the .tsv.gz format -> lets convert it into a csv file to use with our bash script----------
convertTSVtoCSV <- function(input_file, fileprefix){
  # Define the output file name
  output_file <- paste0(fileprefix, ".csv")
  
  # Define the output file path (if you want to specify a directory)
  output_file_path <- paste0("/Volumes/Expansion/data1/Day-25+55_whitelists-SR/", output_file)
  
  # Read the gzipped TSV file
  data <- read_tsv(input_file, show_col_types = FALSE)
  
  # Write the data to a CSV file
  write_csv(data, output_file_path)
  
  return(output_file_path)
}

setwd("/Volumes/Expansion/data1/Day-25+55_whitelists-SR")
sr.D25.tsv.path <- "/Volumes/Expansion/CELLRANGER_counts/sr_day25_cortdiff/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
sr.D55.tsv.path <- "/Volumes/Expansion/CELLRANGER_counts/sr_day55_cortdiff/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

sr.D25.csv.path <- convertTSVtoCSV(sr.D25.tsv.path, "sr-cortdiff-day25-whitelist")
sr.D55.csv.path <- convertTSVtoCSV(sr.D55.tsv.path, "sr-cortdiff-day55-whitelist")


# Call the compare barcode bash script I wrote using system2()----------------
setwd("/Volumes/Expansion/data1/barcodeComparison")
system2("/Volumes/Expansion/Scripts/compareWhitelists.sh", args = c(lr.D25.barcodes.path, sr.D25.csv.path))
system2("/Volumes/Expansion/Scripts/compareWhitelists.sh", args = c(lr.D55.barcodes.path, sr.D55.csv.path))

compareBarcodesfromSeurat <- function(short.read.seurat.obj, long.read.seurat.obj,
                                      sampleTimepoint){
  filtered.sr.D25.barcode.df <- as.data.frame(rownames(short.read.seurat.obj@meta.data))
  filtered.lr.D25.barcode.df <- as.data.frame(rownames(long.read.seurat.obj@meta.data))
  View(filtered.lr.D25.barcode.df)
  filtered.lr.D25.barcode.df[, 1] <- gsub('\\.', '-', filtered.lr.D25.barcode.df[, 1])
  
  setwd('/Volumes/Expansion/data1/barcodeComparison')
  longRead.csv.filename <- paste0("longReadBarcodes", sampleTimepoint, ".csv")
  shortRead.csv.filename <- paste0("shortReadBarcodes", sampleTimepoint, ".csv")
  
  write.table(filtered.lr.D25.barcode.df, longRead.csv.filename, sep = ",", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(filtered.sr.D25.barcode.df, shortRead.csv.filename, sep = ",", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  system2("/Volumes/Expansion/Scripts/compareWhitelists.sh", args = c(longRead.csv.filename, shortRead.csv.filename))
  
  setwd("/Volumes/Expansion/Scripts/scRNAseq-Analysis-pipeline")
}

compareBarcodesfromSeurat(short.read.day25, long.read.day25,
                          sampleTimepoint = "D25")
compareBarcodesfromSeurat(short.read.day55, long.read.day55,
                          sampleTimepoint = "D55")


# Create a stacked bar chart to visualize the barcode overlap between SR and LR scRNAseq (both samples pre and post Seurat QC)--------------

# Function to create a data frame for ggplot2 from our bash script
createBarcodeDataframe <- function(unique.barcode.counts, total.barcode.counts, 
                                   sampleLabel, beforeQC = TRUE){
  data <- data.frame(
    File = c("BLAZE\n(Long Reads)", "Cellranger\n(Short Reads)"),
    Unique_Barcodes = unique.barcode.counts,
    Total_Barcodes = total.barcode.counts
  )
  
  # Calculate overlapping barcodes (Total - Unique)
  data$Overlapping_Barcodes <- data$Total_Barcodes - data$Unique_Barcodes
  
  # Reshape data for ggplot2
  data_melted <- reshape2::melt(data, id.vars = "File")
  data_melted <- subset(data_melted, variable != "Total_Barcodes")
  
  data_melted$variable <- as.character(data_melted$variable)
  data_melted[1, "variable"] <- "Unique BLAZE barcodes"
  data_melted[2, "variable"] <- "Unique Cellranger barcodes"
  data_melted[3, "variable"] <- "Overlapping Barcodes"
  if(beforeQC){
    data_melted$Sample <- paste(sampleLabel, "(Before QC)")
  }else{
    data_melted$Sample <- paste(sampleLabel, '(After QC)')
  }
  data_melted <- data_melted[-4, ]
  data_melted <- data_melted[ , -1]
  data_melted$variable <- as.factor(data_melted$variable)
  View(data_melted)
  
  
  data_melted
}


# Create a dataframe with information from day 25 sample barcode comparison bash script outputs
D25.barcode.comparison.df <- createBarcodeDataframe(unique.barcode.counts = c(14, 158),
                                                    total.barcode.counts = c(811, 955),
                                                    sampleLabel = "Day 25")
# Create a dataframe with information for day 55 sample
D55.barcode.comparison.df <- createBarcodeDataframe(unique.barcode.counts = c(0, 169),
                                                    total.barcode.counts = c(336, 505),
                                                    sampleLabel = "Day 55")

# Create a dataframe with information for day 25 sample after all the seurat QC that was performed
QCed.D25.barcode.comparison.df <- createBarcodeDataframe(unique.barcode.counts = c(13, 71),
                                                         total.barcode.counts = c(686, 744),
                                                         sampleLabel = "Day 25",
                                                         beforeQC = FALSE)
# Create a dataframe with information for day 55 sample after all the seurat QC that was performed
QCed.D55.barcode.comparison.df <- createBarcodeDataframe(unique.barcode.counts = c(9, 39),
                                                         total.barcode.counts = c(274, 304),
                                                         sampleLabel = "Day 55",
                                                         beforeQC = FALSE)



combined.barcode.df <- rbind(D25.barcode.comparison.df, QCed.D25.barcode.comparison.df,
                             D55.barcode.comparison.df, QCed.D55.barcode.comparison.df)
# Convert the 'Sample' column to a factor and specify the order that is shown in the plot
combined.barcode.df$Sample <- factor(combined.barcode.df$Sample, 
                                     levels = c("Day 25 (Before QC)",
                                                "Day 25 (After QC)",
                                                "Day 55 (Before QC)", 
                                                "Day 55 (After QC)"))
#View(combined.barcode.df)

# Create the plot
ggplot(combined.barcode.df, aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Blues") +
  ylab("Number of Barcodes") +
  xlab("Sample") +
  ggtitle("Cell barcode overlap between Cellranger and BLAZE") +
  theme_minimal()


#display.brewer.all(colorblindFriendly = TRUE)      # Show all color palettes
