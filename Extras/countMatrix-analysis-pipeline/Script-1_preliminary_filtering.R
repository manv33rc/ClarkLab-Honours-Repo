# Load libraries
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)

# Specify the prefix of the output filtered matrix csv file
output_prefix = "q20_Day25"


## Read in the transcript count csv file----
count_matrix <- read.csv("/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_1000expCells/transcript_count.csv",
                         row.names = 1)

# Remove gene_id column from count matrix
count_matrix <- subset(count_matrix, select = -gene_id)
#transcript_gene_key <- read.csv("/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_1000expCells/transcript_count.csv")
#transcript_gene_key <- transcript_gene_key[,c("transcript_id", "gene_id")]

sum(count_matrix) # find total feature counts across all cells
# 1261683

## For each transcript, calculate the percentage of cells that express it
transcript.percent.expression <- rowMeans(count_matrix > 0) * 100
transcript.percent.expression
# Convert the calculated percentages into a data frame
df <- as.data.frame(transcript.percent.expression)

# Create a boxplot to visualize the distribution of transcript expression percentages across cells for each gene
before.filter.box <- ggplot(df, aes(y = "", x = transcript.percent.expression)) +
  geom_boxplot() +
  labs(y = "", x = "Isoform expression percentage (%)", title = "BEFORE removing isoforms \nexpressed in =<1% of cells") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_discrete(breaks = NULL)
#Generate a table of 5 point summary stastics for transcript expression percentages across cells for each gene
summary_table.before <- as.data.frame(summary(df)) %>% subset(., select=Freq) 
colnames(summary_table.before)[1] <- "Isoform expression frequency \nacross cells (%) \nBEFORE Filtering (%)"


## Filter out transcripts that are expressed in <=1% of all cells (to avoid potential confounding false positives from FLAMES)----
transcript.filter <- names(transcript.percent.expression[transcript.percent.expression>1])  #select isoforms expressed in at >= 1% of cells
count_matrix.subset <- count_matrix[transcript.filter,] # create subset of our count matrix with only the aforementioned isoforms
num_removed_trans = length(names(transcript.percent.expression)) - length(transcript.filter)

table1 <- data.frame(
  "Number of Isoforms Removed" = num_removed_trans,
  "Total Reads Before Filtering" = sum(count_matrix),
  "Total Reads After Filtering" = sum(count_matrix.subset),
  "Number of Raw Read Counts Removed" = sum(count_matrix) - sum(count_matrix.subset)
) %>% 
  t(.)

# plot filtered transcript count matrix to view distribution of transcript expression across cells
filtered.trans.percent.expression <- rowMeans(count_matrix.subset > 0 ) * 100   
df2 <- as.data.frame(filtered.trans.percent.expression)
after.filter.box <- ggplot(df2, aes(y = "", x = filtered.trans.percent.expression)) +
  geom_boxplot() +
  labs(y = "", x = "Expression Percentage (%)", title = "AFTER removal") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_discrete(breaks = NULL)

# Generate a table of 5 point summary stastics for transcript expression percentages 
# across cells for each gene
summary_table.after <- as.data.frame(summary(df2)) %>% subset(., select=Freq) 
colnames(summary_table.after)[1] <- "Isoform expression frequency \nAFTER filtering (%)"

summary_table.concat <- cbind(summary_table.before, summary_table.after)

#arrange all filtering summary statistics into a grid
summary_statistics.figs <- grid.arrange(before.filter.box, after.filter.box, 
             tableGrob(table1), tableGrob(summary_table.concat), nrow = 2)

### Initialize a Seurat objects with raw (non-normalized) data (before and after filtering)----
seurat.obj <- CreateSeuratObject(count_matrix)
seurat.obj.filtered <- CreateSeuratObject(count_matrix.subset)
# seurat.obj.filtered
# 13600 features across 811 samples within 1 assay 

# Visualize RNA molecule count and transcript counts (in metadata) to assess the quality of cells after filtering
vlnplot1 <- VlnPlot(seurat.obj, features = c("nCount_RNA")) + labs(title = "Total RNA molecules (reads)\nper cell BEFORE filtering") + 
  theme(legend.position = "none")
vlnplot2 <- VlnPlot(seurat.obj, features = c("nFeature_RNA")) + labs(title = "Total unique isoforms\nper cell BEFORE filtering") + 
  theme(legend.position = "none")
association.plot1 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + labs(title = "Correlation plot of RNA molecules and distinct transcripts per cell BEFORE filtering") + 
  theme(legend.position = "none")

vlnplot3 <- VlnPlot(seurat.obj.filtered, features = c("nCount_RNA")) + labs(title = "Total RNA molecules (reads)\nper cell AFTER filtering") + 
  theme(legend.position = "none")
vlnplot4 <- VlnPlot(seurat.obj.filtered, features = c("nFeature_RNA")) + labs(title = "Total unique isoforms\nper cell AFTER filtering") + 
  theme(legend.position = "none")
association.plot2 <- FeatureScatter(seurat.obj.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm") + labs(title = "Correlation plot of RNA molecules and distinct transcripts per cell AFTER filtering") + 
  theme(legend.position = "none")

# Violin plot comparisons
vlnplot_comparisons <- grid.arrange(vlnplot1, vlnplot3, vlnplot2, vlnplot4, nrow = 2)
association_comparisions <- grid.arrange(association.plot1, association.plot2, nrow = 2)

# Export the filtered matrix so that it can be used in downstream analyses
filename <- paste("filt_", output_prefix, "_transcript_count.csv", sep = "")
write.csv(count_matrix.subset, filename)

# Export all figures generated in script as a pdf
concatenated.layout <- grid.arrange(summary_statistics.figs, vlnplot_comparisons, association_comparisions,
                                    top=textGrob("CHECK THAT DISTRIBUTION OF READS AND ISOFORMS\nARE GOOD AFTER FILTERING FEATURES\n(to reduce false positives from FLAMES isoform quantification)\n"))
pdfFilename <- paste(output_prefix, "-script1.pdf", sep="")
ggsave(pdfFilename, concatenated.layout, width = 11, height = 25)

# if all cells seem to be good quality based on the preliminary feature filter -
# move on to more QC in next script

