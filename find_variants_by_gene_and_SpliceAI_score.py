## Script to do initial PanelApp filtering to find variants in genes of interest with SpliceAI scores >=0.2
## Input files are specified on the command line and are as follows:
## --samples: Tab separated list of "ID	vcf location" - should include full path to VCFs (see example_sample_file.txt)
## example running command:
## python find_variants_by_gene_and_SpliceAI_score.py --samples example_sample_file.txt
## NB - panel_file = open("/path/to/panel_file.tsv")  needs to be changed to your actual panel file
## NB - SAI_snvs = gzip.open("/path/to/spliceai_scores.masked.snv.hg38.vcf.gz") needs to be changed to your SpliceAI file location
## NB - SAI_indels = gzip.open("/path/to/spliceai_scores.masked.indel.hg38.vcf.gz") needs to be changed to your SpliceAI file location

import gzip
import os
import argparse

def get_options():
	parser = argparse.ArgumentParser(description="###")
	parser.add_argument("--samples", required=True, help="# samples file to process")
	args = parser.parse_args()
	return args
args = get_options()

## Set up input and output files:
infile_samples = open(args.samples)
Outfile = ''.join((args.samples, "_variants_out_SpliceAI.txt")) ## output will be named after the samples file specified with _variants_out_SpliceAI.txt appended
outfile = open(Outfile, 'w')

def get_format(format_string):
    """ figure out the format of the sample columns, map data to position
    Args:
        format_string: text from format column of VCF, colon separated string
    """
    format = format_string.split(":")
    format = dict(zip(format, range(0, len(format))))
    return(format)


## Go through the SpliceAI files and store variants in genes of interest with SpliceAI scores >= 0.2 (can be customised)
gene_dict = {} 
gene_panel = open("/path/to/panel_file.tsv") ## change this to match the panel you want to use

for line in gene_panel:
	line = line.strip()
	words = line.split('\t')
	if len(words) <14: continue ## skip lines missing info (added to address an error where some files had incomplete lines)
	if words[1] != "gene": continue ## skip non-gene entries 
	if words[7].startswith("MONOALLELIC") or words[7].startswith("BOTH") or words[7].startswith("BIALLELIC"): ## this can be modified depending on the type of genes you want to include
		gene_dict[words[2]] = 0

SAI_snvs = gzip.open("/path/to/spliceai_scores.masked.snv.hg38.vcf.gz") ## change this to reflect your SpliceAI file location
SAI_indels = gzip.open("/path/to/spliceai_scores.masked.indel.hg38.vcf.gz") ## change this to reflect your SpliceAI file location

SAI_dict = {}

## SpliceAI SNVs file
for line in SAI_snvs:
	line = line.strip()
	words = line.split('\t')
	if line.startswith('#'): continue ## skip headers
	info = words[7].split('|') ## split the info for parsing: 2-5 are scores, 6-9 are locations
	gene = info[1]
	if gene not in gene_dict: continue ## skip entries that aren't in genes of interest
	max = 0.00 ## set baseline maximum to 0 to compare scores against
	for i in info[2:6]: ## for each of the scores
		if float(i) > float(max): ## see if the score is higher than the current max
			max = i ## if current score is higher than current max, store the score as max
	if float(max) >= 0.2: ## if the maximum score is greater than 0.2 we'll output the variant to a temporary SpliceAI subset which can be deleted after
		variant = ''.join(("chr", words[0], "-", words[1], "-", words[3], "-", words[4]))
		SAI_dict[variant] = line

## SpliceAI indels file
for line in SAI_indels:
	line = line.strip()
	words = line.split('\t')
	if line.startswith('#'): continue ## skip headers
	info = words[7].split('|') ## split the info for parsing: 2-5 are scores, 6-9 are locations
	gene = info[1]
	if gene not in gene_dict: continue ## skip entries that aren't in genes of interest
	max = 0.00 ## set baseline maximum to 0 to compare scores against
	for i in info[2:6]: ## for each of the scores
		if float(i) > float(max): ## see if the score is higher than the current max
			max = i ## if current score is higher than current max, store the score as max
	if float(max) >= 0.2: ## if the maximum score is greater than 0.2 we'll output the variant to a temporary SpliceAI subset which can be deleted after
		variant = ''.join(("chr", words[0], "-", words[1], "-", words[3], "-", words[4]))
		SAI_dict[variant] = line

## Go through VCFs and see if any probands have variants which were stored from SpliceAI files
for line in infile_samples:
	line = line.strip()
	words = line.split('\t')
	if line.startswith('Participant'): continue ## skip header if present
	ID = words[0]
	vcf_file_loc = words[1]
	if os.path.exists(vcf_file_loc): ## Check the VCF exists before trying to open it
		vcf_file = gzip.open(vcf_file_loc) ## If VCFs aren't gzipped, remove "gzip."
		for Line in vcf_file:
			Line = Line.strip()
			Words = Line.split('\t')
			if Line.startswith('#'): continue ## Skip vcf headers
			if Words[6] != "PASS": continue ## Skip anything that doesn't have a PASS in the filter column
			if Words[9].startswith('0/0'): continue ## Skip anything where the proband doesn't actually have a variant here
			variant = '-'.join((Words[0], Words[1], Words[3], Words[4]))
			if variant not in SAI_dict: continue
			format = get_format(Words[8]) ## allows it to pull out the DP from the genotype field
			get_DP = Words[9].split(':')
			DP = get_DP[format["DP"]] 
			if int(DP) < 6: continue ## check the depth is at least 5 reads (can be changed)
			outline = ''.join((ID, '\t', Line, '\t', SAI_dict[variant], '\n'))
			outfile.write(outline)
	else:
		print "File ", vcf_file_loc, "not found"
outfile.close()



