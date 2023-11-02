#!/bin/bash

# Script that compares the barcodes listed in two BLAZE-generated whitelist.csv files and returns barcodes that don't match, to run the script
# navigate to the directory that contains the compareWhitlists.sh file
# then execute it using ./compareWhitelists.sh [~/file1.csv] [~/file2.csv] in the terminal

# The script will generate two output files :
# 1) a text file that contains the paths that correspond to each file identifier used in the csv file (File 1 and File 2), 
# and the number non-matching barcodes contained within those files
# 2) a csv file that lists the sequences of the non-matching barcodes with their respective file identifiers

file1="$1"
file2="$2"

# Assign the name of the txt containing the non matching barcode count for each file
# and the csv file listing the actual non matching barcodes with their file identifier
txtFileOutput="Non Matching Barcodes Counts.txt"
csvFileOutput="List of non-matching barcodes.csv"

# Write the header of the csv file
echo "Barcode,File of Origin" > "$csvFileOutput"

# Check if both file paths are provided
if [ -z "$file1" ] || [ -z "$file2" ]; then
  echo "Please provide two BLAZE whitelist.csv file paths as parameters."
  exit 1
fi

# Get the number of lines in each file
file1_lines=$(wc -l < "$file1")
file2_lines=$(wc -l < "$file2")

# Compare the lines of the two CSV files and write non-matching lines to the output file
comm_output=$(comm -23 <(sort "$file1") <(sort "$file2"))
non_matching_count_file1=$(comm -23 <(sort "$file1") <(sort "$file2") | wc -l)
non_matching_count_file2=$(comm -13 <(sort "$file1") <(sort "$file2") | wc -l)

# Write the file paths, line counts, and non-matching line counts to the output file
echo "File 1 ($file1_lines barcodes): $file1" > "$txtFileOutput"
echo "File 2 ($file2_lines barcodes): $file2" >> "$txtFileOutput"
echo "Non-matching barcodes found in File 1: $non_matching_count_file1" >> "$txtFileOutput"
echo "Non-matching barcodes found in File 2: $non_matching_count_file2" >> "$txtFileOutput"

echo
echo "Non-matching barcode counts for each file + file identifiers are saved in $txtFileOutput"

# Append the non-matching lines with the corresponding file identifier to the csv output file:
while read -r line; do
  echo "$line,File 1" >> "$csvFileOutput"
done < <(comm -23 <(sort "$file1") <(sort "$file2"))

while read -r line; do
  echo "$line,File 2" >> "$csvFileOutput"
done < <(comm -13 <(sort "$file1") <(sort "$file2"))

echo "List of non-matching barcodes are saved in $csvFileOutput"
echo
