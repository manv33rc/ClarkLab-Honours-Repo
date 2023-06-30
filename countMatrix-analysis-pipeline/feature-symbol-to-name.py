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
    elif featureConversionType == 'trans.to.gene':
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

def FeatureSymboltoName(feature_id, ConvType, gtfFilePath):
    valid_ConversionTypes = ['transcript', 'gene', 'trans.to.gene']
    if ConvType not in valid_ConversionTypes:
        print("Error: Invalid Conversion Type \nRequires: 'transcript', 'gene' or 'trans.to.gene'")
    
    # Build a hashmap of feature symbol and name key-value pairs
    feature_map = build_feature_map(ConvType, gtfFilePath)
    # Retrieve the gene name based on gene ID
    feature_name = feature_map.get(feature_id)
    
    if feature_name:
        if ConvType == 'transcript':
            print(f"Transcript symbol: {feature_id}, Transcript Name: {feature_name}")
        elif ConvType == 'gene':
            print(f"Gene symbol: {feature_id}, Gene Name: {feature_name}")
        elif ConvType == 'trans.to.gene':
            print(f"Transcript symbol: {feature_id}, Gene Name: {feature_name}")
    else:
        print(f"No matching names found for the provided feature symbol: {feature_id}")
    return feature_name

# TEST FUNCTIONS
gtf_file = '/Users/manveerchuahan/SCRIPTS/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
gene_id = 'ENSG00000223972.6'   # Gene Name: DDX11L1
transcript_id = 'ENST00000641515.2'  # transcript name : OR4F5-201
#gene_id = 'ENSG00000186092.7'   # Gene Name: OR4F5

# Retrieve desired name based on feature symbol and conversion type
FeatureSymboltoName(gene_id, 'gene', gtf_file)
FeatureSymboltoName(transcript_id, 'transcript', gtf_file)
FeatureSymboltoName(transcript_id, 'trans.to.gene', gtf_file)