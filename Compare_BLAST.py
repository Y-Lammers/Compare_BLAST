#!/usr/bin/python

# Usage: Compare_Sample.py [Output] [Abundance] [Taxon Level] [Sample list]

# Import relevant modules
import sys, copy, os


def GetTaxonomy():

	# Parse through the sample list and obtain the taxonomies present

	# The empty taxon_dic where all taxonomies will be stored
	taxon_dic = {}

	# For each sample
	for sample in sys.argv[4]:
		# for each blast hit in the sample
		for hit in open(sample):
			# Skip if it is the header line
			if 'Query' in hit: continue

			# split into columns
			hit = hit.strip().split('\t')

			# Check if abundance data is present
			# and extract the correct column containing
			# the taxonomy data			
			if sys.argv[2] == '-A':
				taxonomy = ' / '.join(hit[11].split(' / ')[:sys.argv[3]])
			else:
				taxonomy = ' / '.join(hit[10].split(' / ')[:sys.argv[3]])

			# if the taxonomy is not in the dictionary
			# add it with size 0
			if taxonomy not in taxon_dic:
				taxon_dic[taxonomy] = 0

	# return the taxonomy dictionary
	return taxon_dic


def CollectSamples(taxon_dic):
	
	# Get the taxonomy hit number per taxonomy / per sample

	# Create an empty dictionary that will contain the
	# samples + taxonomy numbers
	sample_dic = {}

	# for each sample
	for sample in sys.argv[4]:
		# strip the directories + extension from the sample
		c_sample = os.path.splitext(os.path.basename(sample))[0]

		# create a new subdictionary for the sample in
		# the sample_dic containing the taxon_dic
		sample_dic[c_sample] = copy.deepcopy(taxon_dic)

		# go through the BLAST hits in the sample
		for hit in open(sample):
			# skip if it is the header line
			if 'Query' in hit: continue

			# split into columns
			hit = hit.strip().split('\t')

			# If abundance is present: add it for each
			# BLAST hit, if not increase by 1
			if sys.argv[2] == '-A':
				sample_dic[c_sample][' / '.join(hit[11].split(' / ')[:sys.argv[3]])] += int(hit[1])
			else:
				sample_dic[c_sample][' / '.join(hit[10].split(' / ')[:sys.argv[3]])] += 1

	# return the sample dictionary
	return sample_dic


def PrintResults(sample_dic, taxon_dic):

	# Write the taxonomies + the abundance numbers per sample

	# the get samples and taxonomies present
	# and order these alphabetical
	sample_keys, taxon_keys = sorted(sample_dic.keys()), sorted(taxon_dic.keys())

	# open the output file
	output_file = open(sys.argv[1], 'w')

	# write the output file header
	output_file.write('Taxonomy\t{0}\n'.format('\t'.join(sample_keys)))

	# for taxonomy in the taxon_key list, get the abundance for each
	# of the sample files and write this to the output file
	for taxonomy in taxon_keys:
		output_file.write('{0}\t{1}\n'.format(taxonomy,
		'\t'.join([str(sample_dic[sample][taxonomy]) for sample in sample_keys])))


# Prepare the input parameters: Condense the samples to a 
# single list and convert the taxon level to an integer
sys.argv[4] = [sample for sample in sys.argv[4:]]
sys.argv[3] = int(sys.argv[3])

# Get the taxonomies in the samples
taxon_dic = GetTaxonomy()

# Get the samples and taxonomy abundances
sample_dic = CollectSamples(taxon_dic)

# Print the results
PrintResults(sample_dic, taxon_dic)
