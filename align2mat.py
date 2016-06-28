# -*- coding: utf-8 -*-

################################### Imports ###################################
import numpy as np
import sys, pickle


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

# Correspondances between nucleotides and character
nt_to_char = {"A" : 0, "T" : 1, "C" : 2, "G" : 3}

# Import list of alignment file to process
if len(sys.argv) >= 2 : # List of files file name as an argument
	LIST_FILE = sys.argv[1].split()[0]
	if len(sys.argv) >= 3 : # Matrices type specified
		# 0 : matrices contain 0, 1, 2 or 3 for A, T, C and G
		# 1 : matrices contain 0 or 1 for identical to the reference sequence
		# or mutated from the reference
		MAT_TYPE = int(sys.argv[2])
		if len(sys.argv) >= 4 : # Coordinates type specified
		# 0 : keep coordinates of the reference allele
		# N : use coordinates of GRChN (for N positive)
			COORD_SYST = int(sys.argv[3])
			if len(sys.argv) >= 5 : # Additional checks specified
				ADD_CHECKS = bool(int(sys.argv[4]))
			else :
				ADD_CHECKS = 0
		else :
			COORD_SYST = 0
	else :
		MAT_TYPE = 0
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
	# Reference sequence
	name_of_ref_seq_in_list = []

	# Parse file
	line_index = 0
	on_block = False # Are we reading a block ?
	
	# Used for additional checks
	ref_full_seq = ''

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
				# If we are parsing the reference allele and we want final positions
				# using GRCh coordinates with additional checks
				if ADD_CHECKS and (COORD_SYST != 0) and (len(name_of_ref_seq_in_list) == 0 or allele_name in name_of_ref_seq_in_list):
					# We need to keep the position of pipes
					ref_full_seq  += ''.join(splitted_line[1:])
				allele_seq = ''.join(splitted_line[1:]).replace('|', '')
				# If allele already in dictionary, add seq
				if allele_name in seq_dict.keys():
					seq_dict[allele_name] += allele_seq
				# Otherwise, create entry
				else :
					seq_dict[allele_name] = allele_seq
					# We then check if it is the first sequence of the file,
					# hence the reference sequence
					if len(name_of_ref_seq_in_list) == 0:
						name_of_ref_seq_in_list.append(allele_name)
			line_index += 1

	print "Reference sequence for "+ file +" is "+name_of_ref_seq_in_list[0]

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
	
	# Create alignment matrix
	align_matrix = np.empty([len(seq_dict.values()), len(pos_list)])
	# Order in which sequences are reported in the matrix
	ordered_seq_list = []
	
	# Start at first line of the matrix
	i = 0
	if MAT_TYPE == 1 : # Fill matrix with 0 or 1
		for allele_name, allele_seq in seq_dict.iteritems():
			ordered_seq_list.append(allele_name)
			if allele_name != name_of_ref_seq_in_list[0]:
				for j in xrange(len(pos_list)):
					if allele_seq[pos_list[j]-1] == "-":
						align_matrix[i][j] = 0
					else:
						align_matrix[i][j] = 1
			else:
				align_matrix[i] = [0]*len(pos_list)
			i += 1

	elif MAT_TYPE == 0 : # Fill matrix with 0, 1, 2 or 3
		for allele_name, allele_seq in seq_dict.iteritems():
			ordered_seq_list.append(allele_name)
			if allele_name != name_of_ref_seq_in_list[0]:
				for j in xrange(len(pos_list)):
					if allele_seq[pos_list[j]-1] == "-":
						align_matrix[i][j] = nt_to_char[ seq_dict[name_of_ref_seq_in_list[0]][pos_list[j]-1] ]
					else:
						align_matrix[i][j] = nt_to_char[ allele_seq[pos_list[j]-1] ]
			else:
				align_matrix[i] = [nt_to_char[allele_seq[pos_list[k]-1]] for k in xrange(len(pos_list))]
			i += 1
	
	else :
		print("Unknown matrix type.")

	# The positions we kept are refering to positions in the alignment and not
	# positions in the reference sequence. The reference sequence can indeed
	# include deletions compared to other sequences. We need to modify the
	# positions list to take this fact into account.

	# List of positions on which a shift occurs.
	# NB : the indices correspond to the number of shifts at this asociated
	# position
	shifts = [0]

	# Find the positions of those shifts
	ref_seq = seq_dict[name_of_ref_seq_in_list[0]]
	for i in xrange(len(seq_dict.values()[0])):
		if ref_seq[i] == ".":
			shifts.append(i+1)

	# i is the index for parsing the shifts list
	i = len(shifts)-1
	corrected_pos_list = []
	# We parse the position list in reverse order
	# In order to have the biggest numbers first
	for j in reversed(xrange(len(pos_list))):
		while pos_list[j] < shifts[i]:
			i -= 1
		corrected_pos_list.append(pos_list[j] - i)
	pos_list = reversed(corrected_pos_list)

	# If necessary, convert positions to GRCh coordinates
	if COORD_SYST:
		# Read biological data from ENSEMBL for this given version of GRCh
		bio_data_file = open(file+".GRCh"+str(COORD_SYST)+".bio", "r")
		bio_data = [int(i.split()[0]) for i in bio_data_file.readlines()]
		
		# Is the gene coding on the direct or reverse strand ?
		direct_strand = bio_data[0]
		bio_data = bio_data[1:]

		base = bio_data[0] # First line correspond to the start position on GRCh
		
		# Index used to parse bio_pos, every positive odd line corresponds
		# to the length of a coding sequence, and every even line to the length
		# of an intron
		i = 1 

		# We want to compare the position on the reference sequence to the cumulated
		# length of exons, so we compute the cumulated sums from the individual lengths
		k = 1
		while 2*k+1 < len(bio_data):
			bio_data[2*k+1] += bio_data[2*(k-1)+1]
			k += 1
		
		corrected_pos_list = []
		# Parse on each kept position
		for pos in pos_list:
			while pos > bio_data[i]: # While we are parsing a sequence after at least 
			# one new intron
				if direct_strand:
					base += bio_data[i+1] # We add the new intron to the shift at each position
				else:
					base -= bio_data[i+1]
				i += 2
			if direct_strand:
				corrected_pos_list.append(pos + base - 1)
			else:
				corrected_pos_list.append(1 + base - pos)
		pos_list = corrected_pos_list

		# ADDITIONAL CHECKS
		if ADD_CHECKS:
			print "Performing additional checks"
			# Check if the pipes in the reference sequence truly corresponds to the introns
			# with respect to data provided on ENSEMBL
			intron_pos_in_file = locations_of_substring(ref_full_seq.replace('.', ''), '|')
			for i in xrange(len(intron_pos_in_file)):
				intron_pos_in_file[i] -= i
			print "Positions of introns in alignment match positions according to ENSEMBL:"
			print intron_pos_in_file == [bio_data[2*k+1] for k in xrange((len(bio_data)-1) / 2)]

	# Write a pickle file including the computed matrix
	pickle.dump(align_matrix, open(file+".mat",'wb'))
	# Write files including the order of the sequences in the matrix
	# And the positions kept, 1-based indices, based on the reference sequence
	order_file = open(file+'.ord', 'w')
	order_file.write(" ".join(ordered_seq_list))
	order_file.close()

	pos_file = open(file+'.pos', 'w')
	pos_file.write(" ".join([str(i) for i in pos_list]))
	pos_file.close()
	file_data.close()
LIST_FILE_DATA.close()