#!/usr/bin/python

# Usage: Compare_Sample.py [Output] [Abundance] [Sample list]

# Import relevant modules
import sys, copy, os


def GetTaxonomy():

	# Parse through the sample list and obtain the taxonomies present

	# The empty taxon_dic where all taxonomies will be stored
	taxon_dic = {'total': 0}

	# For each sample
	for sample in sys.argv[3]:
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
				taxonomy = hit[11]
			else:
				taxonomy = hit[10]

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
	for sample in sys.argv[3]:
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
				sample_dic[c_sample][hit[11]] += int(hit[1])
				sample_dic[c_sample]['total'] += int(hit[1])
			else:
				sample_dic[c_sample][hit[10]] += 1
				sample_dic[c_sample]['total'] += 1

	# return the sample dictionary
	return sample_dic


def PrintResults(sample_dic, taxon_dic):

	# Write the taxonomies + the abundance numbers per sample

	# the get samples and taxonomies present
	# and order these alphabetical
	sample_keys, taxon_keys = sorted(sample_dic.keys()), sorted(taxon_dic.keys())

	# open the output file
	output_file = open(sys.argv[1], 'w')

	# Prepare the taxonomic header
	taxon_head = 'Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies'
	sample_head = '\t'.join(['{0}\tRelative Abundance'.format(sample) for sample in sample_keys])

	# write the output file header
	output_file.write('{0}\t{1}\n'.format(taxon_head, sample_head))

	# for taxonomy in the taxon_key list, get the abundance (both raw and
	# relative abundance) for each of the sample files and write
	# this to the output file
	for taxonomy in taxon_keys:
		if taxonomy == 'total': continue

		# taxonomy data
		tax_format = taxonomy.split(' / ')

		# fill if columns are missing
		while len(tax_format) < 7:
			tax_format.append(' ')

		output_file.write('{0}\t{1}\n'.format('\t'.join(tax_format),'\t'.join(
		["{0}\t{1:.3f}".format(str(sample_dic[sample][taxonomy]),(
		float(sample_dic[sample][taxonomy])/sample_dic[sample]['total'])*100) 
		for sample in sample_keys])))


# Prepare the input parameters: Condense the samples to a 
# single list and convert the taxon level to an integer
sys.argv[3] = [sample for sample in sys.argv[3:]]

# Get the taxonomies in the samples
taxon_dic = GetTaxonomy()

# Get the samples and taxonomy abundances
sample_dic = CollectSamples(taxon_dic)

# Print the results
PrintResults(sample_dic, taxon_dic)
