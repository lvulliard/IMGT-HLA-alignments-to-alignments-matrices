# -*- coding: utf-8 -*-

################################### Imports ###################################
import numpy as np
import sys


################################## Functions ##################################
def locations_of_substring(string, substring):
    """Return a list of locations of a substring."""

    substring_length = len(substring)    
    def recurse(locations_found, start):
        location = string.find(substring, start)
        if location != -1:
            return recurse(locations_found + [location], location+substring_length)
        else:
            return locations_found

    return recurse([], 0)

##################################### Main ####################################

# Import list of alignment file to process
if len(sys.argv) >= 2 : # List of files file name as an argument
	LIST_FILE = int(sys.argv[1])
else :
	LIST_FILE = "alignments_file_list.txt"

# Put alignment files names in an array (names should not include spaces)
LIST_FILE_DATA = open(LIST_FILE)
FILE_NAMES = (name.split()[0] for name in LIST_FILE_DATA.readlines())

# Loop on every file to process
for file in FILE_NAMES:
	print "Processing "+file

	# Read file, skipping header
	file_data = open(file)
	lines = file_data.readlines()[6:]

	# Parsing dictionary
	# Key = Sequence name
	# Value = Sequence
	seq_dict = dict()

	# Parse file
	line_index = 0
	on_block = False # Are we reading a block ?
	while(line_index < len(lines)):
		line = lines[line_index]
		splitted_line = line.split()

		# If we are not already in an alignment block
		if not on_block :
			# If we are on the first line of a block
			if len(splitted_line) > 0 and splitted_line[0] == "cDNA":
				# Go to first allele, either 2 or 3 lines down, depening on file.
				if len(lines[line_index +2].split()) > 1:
					line_index +=2
				else :
					line_index +=3
				on_block = True # We start reading a block
			else:
				# Go to next line
				line_index +=1
		else :
			# If we reach the end of a block
			if len(splitted_line) == 0:
				on_block = False
			else :
				# We are on a line of the alignment
				allele_name = splitted_line[0]
				allele_seq = ''.join(splitted_line[1:]).replace('|', '')
				# If allele already in dictionary, add seq
				if allele_name in seq_dict.keys():
					seq_dict[allele_name] += allele_seq
				# Otherwise, create entry
				else :
					seq_dict[allele_name] = allele_seq
			line_index += 1

	# Data curation
	# We remove positions with unknown nucleotide for a sequence
	# We also remove indels

	# List of positions kept, 1-based indices
	pos_list = range(1, len(seq_dict.values()[0])+1)

	# Parse file
	line_index = 0
	on_block = False # Are we reading a block ?
	for seq in seq_dict.values():
		# If * character in a position kept at this stage
		# It corresponds to an undetermined base in the sequence
		if ''.join([seq[i-1] for i in pos_list]).count('*') > 0 :
			# We remove the corresponding positions from the list
			pos_list = list(set(pos_list) - set([i+1 for i in locations_of_substring(seq, '*')]))
		# We do the same for the . character
		# It corresponds to an indel compared to the reference sequence
		if ''.join([seq[i-1] for i in pos_list]).count('.') > 0 :
			pos_list = list(set(pos_list) - set([i+1 for i in locations_of_substring(seq, '.')]))
