#!/bin/bash
#SBATCH --partition="physical"
#SBATCH --nodes=1
#SBATCH --account="punim0646"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mschauhan@student.unimelb.edu.au
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#time in days-hours:mins:sec
#SBATCH --time=1-00:00:00


# Load conda module then activate the environment that contains minimap2
module load anaconda3/2021.11
source activate minimap2
echo "Conda environment initialized"
# Load module for R
module load foss/2021b r/4.2.0

# Assign file paths
fastqFile="/data/gpfs/projects/punim0646/manveer/BLAZE_paper_fastqs/LR_sc-RNA-seq_GridION_Q20_pass.fastq.gz"
whitelist="/data/gpfs/projects/punim0646/manveer/BLAZE_Q20_15.06/500expectedcells/whitelist.csv"
outputDirectory="/data/gpfs/projects/punim0646/manveer/FLAMES_Q20_GridION_500cells"
RscriptPath="/data/gpfs/projects/punim0646/manveer/Scripts/FLAMES.r"

# Run FLAMES.R
Rscript ${RscriptPath} ${fastqFile} ${whitelist} ${outputDirectory}