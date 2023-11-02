#!/bin/bash
#SBATCH --partition="physical"
#SBATCH --nodes=1
#SBATCH --account="punim0646"
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mschauhan@student.unimelb.edu.au
#SBATCH --cpus-per-task=12
#SBATCH --mem=100000
#time in days-hours:mins:sec
#SBATCH --time=1-00:00:00


# Load conda module then activate the environment that contains minimap2
module load anaconda3/2021.11
source activate minimap2
echo "Conda environment initialized"
# Load module for R
module load foss/2021b r/4.2.0

# Assign file paths
FLAMESpath="/data/gpfs/projects/punim0646/manveer/Scripts/FLAMES.r"

# Run FLAMES.R - Make sure FLAMES.r and config.json are correct
Rscript ${FLAMESpath}