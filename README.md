# Extract_variants
Scripts to extract variants from VCFs when they meet customisable criteria

find_variants_by_gene_and_consequence.py - takes a list of samples/VCFs, panels and genes and extracts variants in those genes that have VEP consequences listed in CQ_dict.

find_variants_by_gene_and_SpliceAI_score.py - takes a list of samples/VCFs, and a panel of genes (to be added in the script) and pulls out variants from the VCFs with SpliceAI scores of 0.2 or greater (you need to point to your copy of the SpliceAI files in the code). 
