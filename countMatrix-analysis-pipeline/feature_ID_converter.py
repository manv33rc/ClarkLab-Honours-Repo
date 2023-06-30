import re

# Function to extract feature information based on conversion type
def extract_feature_info(featureConversionType, line):
    if featureConversionType == 'transcript':
        # Define the regex patterns if matching transcript id to transcript name
        feature_id_pattern = r'transcript_id\s+"([^"]+)"'
        feature_name_pattern = r'transcript_name\s+"([^"]+)"'
    elif featureConversionType == 'gene':
        # Define the regex patterns if matching gene id to gene name
        feature_id_pattern = r'gene_id\s+"([^"]+)"'
        feature_name_pattern = r'gene_name\s+"([^"]+)"'
    elif featureConversionType == 'trans.to.gene':
        # Define the regex patterns if matching transcript id to gene name
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


# Function to build a hashmap of feature symbol and name key-value pairs
def build_feature_map(featureConversionType, gtf_file):
    feature_map = {}  # Initialize the hashmap to store feature ID-name key-value pairs

    with open(gtf_file, 'r') as file:
        for line in file:
            # Extract transcript ID and transcript name from each line
            feature_id, feature_name = extract_feature_info(featureConversionType, line)
            if feature_id and feature_name:
                # Add gene ID and gene name as key-value pair to the hashmap
                feature_map[feature_id] = feature_name

    return feature_map

# Function to retrieve the feature name based on gene ID or transcript ID
def FeatureIDtoSymbol(feature_id, ConvType, gtfFilePath, novelIsoTitle = False):
    # Check if conversion type parameter is in correct format
    valid_ConversionTypes = ['transcript', 'gene', 'trans.to.gene']
    if ConvType not in valid_ConversionTypes:
        print("Error: Invalid Conversion Type \nRequires: 'transcript', 'gene' or 'trans.to.gene'")
    
    # Check if feature_id is a novel isoform
    if '-' in feature_id and novelIsoTitle and ConvType in ['transcript', 'trans.to.gene']:
        novel_isoform_id = feature_id
        feature_id = feature_id.split('-')[0]  # Extract the substring before the first hyphen (the novel isoform's gene)
        # Build hashmap of gene symbol and name pairs and find the novel isoform's corresponding gene name
        feature_map = build_feature_map('gene', gtfFilePath)
        feature_name = feature_map.get(feature_id)
        
        if feature_name:
            fullIsoTitle = f"Novel {feature_name} Isoform \n({novel_isoform_id})"
            print(fullIsoTitle)
            return fullIsoTitle
        else:
            print(f"No matching genes found for novel Isoform ID: {feature_id}")
            return novel_isoform_id

    else:
        # Build a hashmap of feature symbol and name key-value pairs
        feature_map = build_feature_map(ConvType, gtfFilePath)
        # Retrieve the feature name based on feature ID
        feature_name = feature_map.get(feature_id)
    
        if feature_name:
            if ConvType == 'transcript':
                print(f"Transcript ID: {feature_id}, Transcript Symbol: {feature_name}")
            elif ConvType == 'gene':
                print(f"Gene ID: {feature_id}, Gene Symbol: {feature_name}")
            elif ConvType == 'trans.to.gene':
                print(f"Transcript ID: {feature_id}, Gene Symbol: {feature_name}")
            return feature_name
        else:
            print(f"No matching names found for the provided feature ID: {feature_id}")
            return feature_id

# TEST FUNCTIONS
gtf_file = '/Users/manveerchuahan/SCRIPTS/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
gene_id = 'ENSG00000104435.14'   # Gene Name: DDX11L1
transcript_id = 'ENSG00000104435.14-79611117-79665011-1'  # transcript name : OR4F5-201
#gene_id = 'ENSG00000186092.7'   # Gene Name: OR4F5

# Retrieve desired name based on feature symbol and conversion type
FeatureIDtoSymbol(gene_id, 'gene', gtf_file, novelIsoTitle=True)
FeatureIDtoSymbol(transcript_id, 'gene', gtf_file, novelIsoTitle=True)
FeatureIDtoSymbol(transcript_id, 'transcript', gtf_file, novelIsoTitle=True)
FeatureIDtoSymbol(transcript_id, 'trans.to.gene', gtf_file, novelIsoTitle=True)