import re

def extract_feature_info(featureConversionType, line):
    if featureConversionType == 'transcript':
        # Define the regex patterns if matching transcript symbol to transcript name
        feature_id_pattern = r'transcript_id\s+"([^"]+)"'
        feature_name_pattern = r'transcript_name\s+"([^"]+)"'
    elif featureConversionType == 'gene':
        # Define the regex patterns if matching gene symbol to gene name
        feature_id_pattern = r'gene_id\s+"([^"]+)"'
        feature_name_pattern = r'gene_name\s+"([^"]+)"'
    elif featureConversionType == 'trans-to-gene':
        # Define the regex patterns if matching transcript symbol to gene name
        feature_id_pattern = r'transcript_id\s+"([^"]+)"'
        feature_name_pattern = r'gene_name\s+"([^"]+)"'

    # Search for feature ID and feature name within the line
    feature_id_match = re.search(feature_id_pattern, line)
    feature_name_match = re.search(feature_name_pattern, line)

    if feature_id_match and feature_name_match:
        # Extract feature ID and feature name from the match
        feature_id = feature_id_match.group(1)
        feature_name = feature_name_match.group(1)
        return feature_id, feature_name
    else:
        return None, None
    

def build_feature_map(featureConversionType, gtf_file):
    feature_map = {}  # Initialize the hashmap to store gene ID-gene name pairs

    with open(gtf_file, 'r') as file:
        for line in file:
            # Extract transcript ID and transcript name from each line
            feature_id, feature_name = extract_feature_info(featureConversionType, line)
            if feature_id and feature_name:
                # Add gene ID and gene name as key-value pair to the hashmap
                feature_map[feature_id] = feature_name

    return feature_map

def GeneSymboltoName():
    gene_map = 
    return gene_name
# Example usage
gtf_file = '/Users/manveerchuahan/SCRIPTS/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'

# Build the gene map by processing the GTF file
gene_map = build_gene_map(gtf_file)

# Retrieve gene name based on gene ID
#gene_id = 'ENSG00000223972.6'   # Gene Name: DDX11L1
gene_id = 'ENSG00000186092.7'   # Gene Name: OR4F5
gene_name = gene_map.get(gene_id)

if gene_name:
    print(f"Gene ID: {gene_id}, Gene Name: {gene_name}")
else:
    print(f"No matching gene found for ID: {gene_id}")