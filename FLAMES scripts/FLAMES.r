# R script used to generate flames output files from BLAZE generated whitelist.csv file
# USAGE (type in terminal) - Rscript FLAMES.R [/path/to/fastq] [/path/to/whitelist.csv] [path/to/outputDirectory]


# Install the FLAMES package if needed
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("FLAMES")


library("FLAMES")


# Define the variables that are unique to the sample-------------------------
fastqfile <- "/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_1000expCells/matched_reads.fastq.gz"
whitelist <- "/data/gpfs/projects/punim0646/manveer/BLAZE_Q20_15.06/whitelist.csv"
output <- "/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_1000expCells"


# DEFINE GLOBAL VARIABLES HERE (these shouldn't change between samples)-----------
minimap2_path <- "~/.conda/envs/minimap2/bin/minimap2"
# Define reference files
GTF <- "/data/gpfs/projects/punim0646/manveer/gencode.v43.basic.annotation.gtf"
transcriptome <- "/data/gpfs/projects/punim0646/manveer/gencode.v43.transcripts.fa"
# Enter the path for the config file
configPath <- "/data/gpfs/projects/punim0646/manveer/Scripts/config.json"



# Run FLAMES-----------------------------------------------------------------
sce <- sc_long_pipeline(fastq = fastqfile, outdir = output, reference_csv = whitelist, annotation = GTF, 
    genome_fa = transcriptome, match_barcode = FALSE, minimap2_dir = minimap2_path, 
    config_file = configPath)
    
message("Complete")