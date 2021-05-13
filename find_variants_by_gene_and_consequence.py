## Script to do initial PanelApp filtering to find variants in genes of interest with particular consequences (can be specified by altering the CQ_dict dictionary)
## Input files are specified on the command line and are as follows:
## --samples: Tab separated list of "ID	vcf location" - should include full path to VCFs (see example_sample_file.txt)
## --panels: Tab separated list of "ID	panel name" - should include full paths to panels files (see example_panel_file.txt)
## --genes: List of other genes of interest - just a list of gene names, one per line (see example_gene_file.txt)
## example running command:
## python find_variants_by_gene_and_consequence.py --samples example_sample_file.txt --panels example_panel_file.txt --genes example_gene_file.txt

import gzip
import os
import argparse

def get_options():
	parser = argparse.ArgumentParser(description="###")
	parser.add_argument("--samples", required=True, help="# samples file to process")
	parser.add_argument("--panels", required=True, help="# file linking sample to panel(s)")
	parser.add_argument("--genes", required=True, help="# additional gene list to check for all samples")
	args = parser.parse_args()
	return args
args = get_options()

## Set up input and output files:
infile_samples = open(args.samples)
infile_panels = open(args.panels)
infile_genes = open(args.genes)

Outfile = ''.join((args.samples, "_variants_out1.txt")) ## output will be named after the samples file specified with _variants_out1.txt appended
outfile = open(Outfile, 'w')

## Store which panels are relevant for which samples in a dictionary so this file is only parsed once
panel_dict = {}
for line in infile_panels:
	line = line.strip()
	words = line.split('\t')
	if words[0] not in panel_dict: ## If this ID hasn't yet been stored as a key,
		panel_dict[words[0]] = words[1] ## Add the ID to the dictionary with this panel as the value
	else: ## If this ID already has an entry
		panel_dict[words[0]] = panel_dict[words[0]] + ';' + words[1] ## Add this panel to that samples entry separated from previous panels by ;

## Any variants with consequences in this dictionary will be pulled out - consequences can be added or removed as needed
CQ_dict = {"stop_gained": 0, "splice_acceptor": 0, "splice_donor": 0, "frameshift": 0, "missense": 0, "splice_region": 0}

tracking = {} ## setting up a dictionary to store variants that have already been output so we don't get duplicate lines in the output

for line in infile_samples:
	line = line.strip()
	words = line.split('\t')
	if line.startswith('Participant'): continue ## skip header if present
	ID = words[0]
	vcf_file_loc = words[1]
	gene_dict = {} ## This is to store all the relevant gene names for that individual
	if ID not in panel_dict: ## If this ID doesn't have any panels stored, print an error then skip it
		print ID, " - panels unknown"
	if ID not in panel_dict: continue 
	panels = panel_dict[ID].split(';')
	for i in panels:
		panel_file = open(i)
		for iline in panel_file:
			iline = iline.strip()
			iwords = iline.split('\t')
			if len(iwords) < 14: continue ## skip lines missing info (added to address an error where some files had incomplete lines)
			if iwords[1] != "gene": continue ## skip non-gene entries
			if "Expert Review Green" not in iwords[3]: continue ## This skips low confidence genes. This can be commented out if those are wanted.
			if iwords[7].startswith("MONOALLELIC") or iwords[7].startswith("BOTH") or iwords[7].startswith("BIALLELIC"): ## this can be modified depending on the type of genes you want to include
				gene_dict[iwords[2]] = 0
	infile_genes = open(args.genes)
	for lineG in infile_genes: ## go through the additional genes file and add any additional gene names to the gene dictionary
		lineG = lineG.strip()
		wordsG = lineG.split('\t')
		gene_dict[wordsG[0]] = 0
	## Now go through the actual VCFs and start finding variants
	if os.path.exists(vcf_file_loc): ## Check the VCF exists before trying to open it
		vcf_file = gzip.open(vcf_file_loc) ## If VCFs aren't gzipped, remove "gzip."
		for Line in vcf_file:
			Line = Line.strip()
			Words = Line.split('\t')
			if Line.startswith('#'): continue ## Skip vcf headers
			if Words[6] != "PASS": continue ## Skip anything that doesn't have a PASS in the filter column
			if Words[9].startswith('0/0'): continue ## Skip anything where the proband doesn't actually have a variant here
			get_info = Words[7].split(';') ## split up the info field for parsing
			count = 0
			for i in get_info:
				if i.startswith("CSQT="): ## Pull out the bit of the info field that's got variant annotation in it
					split_info = i.split(',')
					count = int(count) +1
			if count == 0: continue ## This skips any lines that don't have VEP variant information
			for i in split_info: ## Go through each part of the split annotation information
				for j in CQ_dict: ## Check if any consequences from our dictionary are in it
					if j in i: ## If this variant's CQ is something we're interested in...
						for k in gene_dict: ## See if it's in a gene we're interested in...
							gene = ''.join(("|",k,"|")) ## Added this in so it matches on complete gene name (before, it would have pulled out anything where a dictionary gene was in another gene's name, e.g. CR1 in dictionary would have pulled out variants in CR1 but also CR1L)
							if gene in i:
								## check the depth is at least 5 (this can be adjusted)
								get_DP = Words[9].split(':')
								DP = get_DP[3] ## NB this will need to be changed if DP is not always in this position
								if int(DP) < 6: continue ## check the depth is at least 5 reads (can be changed)
								## Store the ID and variant to prevent duplicates in the output
								ID_var = '-'.join((words[0], Words[0], Words[1], Words[3], Words[4]))
								if ID_var in tracking: continue ## Skips any entries that have already been stored/output
								tracking[ID_var] = 0
								outline = ''.join((words[0], '\t', Line, '\t', k, '\t', j, '\n')) ## this outputs the ID, the full VCF line, the gene and the consequence
								outfile.write(outline)
	else:
		print "File ", vcf_file_loc, "not found" ## prints an error if the VCF doesn't exist which allows it to carry on processing other samples rather than failing
outfile.close()





