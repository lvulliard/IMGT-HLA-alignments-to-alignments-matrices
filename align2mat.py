#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################### Help ####################################
"""Convert IMGT-HLA alignments to alignment matrices.

Usage:
  align2mat.py [options]
  align2mat.py --help
  align2mat.py --version

Options:
  -l <path>, --list <path>    Path to a list of files to process. [default: alignments_files_list.txt]         
  -c <GRCh>, --coord <GRCh>   Specify output coordinates system. [default: 0]
  -f --filter
  -a --addchecks              Perform additional checks on coordinates change. Only apply when -c is specified.
  -h --help                   Show this screen.
  -v --version                Show version.

"""


################################### Imports ###################################
from docopt import docopt
from collections import defaultdict
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


def position_of_deletion(i, del_pos, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, del_len, original_seq, pos_list):
	try:
		j = pos_list.index(del_pos)
		align_matrix[i][j] = str(align_matrix[i][j])+"D"+str(del_len)
	except ValueError:
		if del_pos in del_pos_list:
			if original_seq[del_pos-1] == ".":
				# We cannot write information to this position since it is deleted
				if direct_strand:
					# Direct strand -> Try to add information to the previous position
					position_of_deletion(i, del_pos-1, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, del_len, original_seq, pos_list)
				else:
					# Reverse strand -> Try to add information to the following position
					position_of_deletion(i, del_pos+1, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, del_len, original_seq, pos_list)
			else:
				indel_pos_to_check[del_pos] += "D"+str(del_len)
		elif del_pos in indel_pos:
			# Deletion is following an insertion, so the deletion must be added before the insertion position
			if direct_strand:
				# Direct strand -> Try to add information to the previous position
				position_of_deletion(i, del_pos-1, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, del_len, original_seq, pos_list)
			else:
				# Reverse strand -> Try to add information to the following position
				position_of_deletion(i, del_pos+1, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, del_len, original_seq, pos_list)


################################### Classes ###################################
class Stretch:
	"""Explicit constructor"""
	def __init__(self, stretch_type, start):
		self.type_insertion = stretch_type
		self.start = start
		self.end = start

	def __str__(self):
		if self.type_insertion:
			return("INS: "+str(self.start)+" - "+str(self.end))
		else:
			return("DEL: "+str(self.start)+" - "+str(self.end))

	def pos(self):
		return(range(self.start, self.end+1))


##################################### Main ####################################

# Correspondances between nucleotides and character
nt_to_char = {"A" : 0, "T" : 1, "C" : 2, "G" : 3}

# Import user input
arguments = docopt(__doc__, version='v0.6')
LIST_FILE = arguments["--list"] # List of files file name
COORD_SYST = int(arguments["--coord"]) # Coordinates type
# 0 : keep coordinates of the reference allele
# N : use coordinates of GRChN
if COORD_SYST != 0:
	ADD_CHECKS = arguments["--addchecks"] # Additional checks needed
else:
	ADD_CHECKS = False
FILTER = arguments["--filter"]

# Put alignment files names in an array (names should not include spaces)
LIST_FILE_DATA = open(LIST_FILE)
FILE_NAMES = [name.split()[0] for name in LIST_FILE_DATA.readlines()]

# If additional checks needed, put reference sequence in an array
LIST_REF = [None]*len(FILE_NAMES)
if ADD_CHECKS:
	LIST_REF_FILE = open(LIST_FILE+".ref"+str(COORD_SYST), "r")
	LIST_REF = [i.split()[0] for i in LIST_REF_FILE.readlines()]
	LIST_REF_FILE.close()

# Loop on every file to process
for file, genome_ref_seq in zip(FILE_NAMES, LIST_REF):
	print "Processing "+file

	# Read file, skipping header
	file_data = open(file)
	lines = file_data.readlines()[6:]

	# Parsing dictionary
	# Key = Sequence name
	# Value = Sequence
	seq_dict = dict()
	# Reference sequence
	alignment_ref_seq = []

	# Parse file
	line_index = 0
	on_block = False # Are we reading a block ?
	
	# Used for additional checks
	aln_ref_full_seq = ''
		
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
				# Skip null alleles
				if not FILTER or allele_name[-1] != "N":
					# If we are parsing the reference allele and we want final positions
					# using GRCh coordinates with additional checks
					if ADD_CHECKS and (len(alignment_ref_seq) == 0 or allele_name in alignment_ref_seq):
						# We need to keep the position of pipes
						aln_ref_full_seq += ''.join(splitted_line[1:])
					
					allele_seq = ''.join(splitted_line[1:]).replace('|', '')
					# If allele already in dictionary, add seq
					if allele_name in seq_dict.keys():
						seq_dict[allele_name] += allele_seq
					# Otherwise, create entry
					else :
						seq_dict[allele_name] = allele_seq
						# We then check if it is the first sequence of the file,
						# hence the reference sequence
						if len(alignment_ref_seq) == 0:
							alignment_ref_seq.append(allele_name)
			line_index += 1

	alignment_ref_seq = alignment_ref_seq[0]
	print "Reference sequence for "+ file +" is "+alignment_ref_seq

	if ADD_CHECKS:
		# Define which allele should be used as reference for defining indels in our matrix
		str_len = len(genome_ref_seq)
		for seq in seq_dict.keys():
			if len(seq) >= str_len and seq[0:str_len] == genome_ref_seq:
				genome_ref_seq = seq
				break

		print "Reference sequence for genome GRCh"+ str(COORD_SYST) +" is "+genome_ref_seq
	
	# Data curation
	# We remove positions with unknown nucleotide for a sequence
	# We also separate indels positions

	# List of positions kept, 1-based indices
	pos_list = range(1, len(seq_dict.values()[0])+1)

	# List of positions with indels
	indel_pos = []

	# List of positions with * symbols
	stars_pos = []

	# Parse file
	line_index = 0
	on_block = False # Are we reading a block ?
	for seq in seq_dict.values():
		# If * character in a position kept at this stage
		# It corresponds to an undetermined base in the sequence
		if ''.join([seq[i-1] for i in pos_list]).count('*') > 0 :
			# We remove the corresponding positions from the lists
			new_stars = [i+1 for i in locations_of_substring(seq, '*')]
			stars_pos = list(set(stars_pos + new_stars))
			pos_list = list(set(pos_list) - set(new_stars))
			indel_pos = list(set(indel_pos) - set(stars_pos))
		# We also check the . character
		# It corresponds to an indel compared to the reference sequence
		if ''.join([seq[i-1] for i in pos_list]).count('.') > 0 :
			new_indels = [i+1 for i in locations_of_substring(seq, '.')]
			indel_pos += new_indels
			indel_pos = list(set(indel_pos) - set(stars_pos))
			pos_list = list(set(pos_list) - set(new_indels))
	
	# Create alignment matrix
	align_matrix = np.empty([len(seq_dict.values()), len(pos_list)], dtype=object)
	# Order in which sequences are reported in the matrix
	ordered_seq_list = []

	# Start at first line of the matrix
	i = 0
	for allele_name, allele_seq in seq_dict.iteritems():
		ordered_seq_list.append(allele_name)
		if allele_name != alignment_ref_seq:
			for j in xrange(len(pos_list)):
				if allele_seq[pos_list[j]-1] == "-":
					align_matrix[i][j] = nt_to_char[ seq_dict[alignment_ref_seq][pos_list[j]-1] ]
				else:
					align_matrix[i][j] = nt_to_char[ allele_seq[pos_list[j]-1] ]
		else:
			align_matrix[i] = [nt_to_char[allele_seq[pos_list[k]-1]] for k in xrange(len(pos_list))]
		i += 1

	# Encode indels

	# Compute stretches positions
	indel_pos = sorted(indel_pos)
	insertions_list = []
	deletions_list = []
	ref_seq = seq_dict[alignment_ref_seq]
	ref_seq_indels = [ref_seq[i-1] for i in indel_pos]
	stretches_list = []
	
	# Start the first stretch on first position
	if ref_seq_indels[0] == ".":
		# We are looking at an insertion compared to the reference
		current_stretch = Stretch(True, indel_pos[0])
	else:
		# We are looking at a deletion compared to the reference
		current_stretch = Stretch(False, indel_pos[0])
	stretches_list.append(current_stretch)

	# Check every indel position
	for i in xrange(1, len(indel_pos)):
		if ref_seq_indels[i] == ".":
			# We are looking at an insertion compared to the reference
			if ref_seq[indel_pos[i]-2] == ".":
				# Previous position was already an insertion
				# We are still in the same stretch, so we shift its end
				current_stretch.end += 1
			else:
				# Previous position was not an insertion
				# We are starting a stretch
				current_stretch = Stretch(True, indel_pos[i])
				stretches_list.append(current_stretch)
		else:
			# We are looking at a deletion compared to the reference
			if indel_pos[i-1] == indel_pos[i]-1 and ref_seq_indels[i-1] != ".":
				# Previous position was already a deletion
				# We are still in the same stretch, so we shift its end
				current_stretch.end += 1
			else:
				# Previous position was not a deletion
				# We are starting a stretch
				current_stretch = Stretch(False, indel_pos[i])
				stretches_list.append(current_stretch)


	# We now add the deletion positions to a matrix
	del_pos_list = []
	for stretch in stretches_list:
		if not stretch.type_insertion:
			del_pos_list += stretch.pos()
	del_align_matrix = np.zeros([len(seq_dict.values()), len(del_pos_list)], dtype=object)
	
	# This matrix is filled the next time we iterate over the sequences
	# And then concatenated to the alignment matrix

	# To know where to store the information on indels
	# We need to know which strand of DNA the gene is on
	if COORD_SYST:
		# Read biological data from ENSEMBL for this given version of GRCh
		bio_data_file = open(file+".GRCh"+str(COORD_SYST)+".bio", "r")
		bio_data = [int(i.split()[0]) for i in bio_data_file.readlines()]
		
		# Is the gene coding on the direct or reverse strand ?
		direct_strand = bio_data[0]
		bio_data = bio_data[1:]
	else:
		# By default we assume the gene is on the direct strand
		# So we encode the indels on the base just before the stretch
		direct_strand = 1

	
	for allele_name, allele_seq in seq_dict.iteritems():
		# Add the information on indels on suitable positions
		
		i = ordered_seq_list.index(allele_name)
		# Store indel position if not found in the kept positions list
		# In case it belongs to the deletions positions
		indel_pos_to_check = defaultdict(str)
		for stretch in stretches_list:
			if stretch.type_insertion:
				ins_seq = (''.join([allele_seq[m-1] for m in stretch.pos()])).replace('.','')
				len_ins_seq = len(ins_seq)
				if len_ins_seq > 0:
					# There is an insertion on this sequence
					if direct_strand:
						# Add information before the stretch
						ins_pos = stretch.start - 1
					else:
						# Add information after the stretch on alignement position, so before in genomics position
						ins_pos = stretch.end + 1
					try:
						j = pos_list.index(ins_pos)
						align_matrix[i][j] = str(align_matrix[i][j])+"I"+str(len_ins_seq)
					except ValueError:
						if ins_pos in del_pos_list:
							indel_pos_to_check[ins_pos] += "I"+str(len_ins_seq)
			else:
				# Fill deletion matrix
				for j in stretch.pos():
					seq_base = allele_seq[j-1]
					if seq_base in ["-", "."] or allele_name == alignment_ref_seq:
						seq_base = nt_to_char[ seq_dict[alignment_ref_seq][j-1] ]
					else:
						seq_base = nt_to_char[ seq_base ]
					del_align_matrix[i][del_pos_list.index(j)] = seq_base

				# Encode deletion
				del_seq = ''.join([allele_seq[m-1] for m in stretch.pos()])
				nb_kept = len(del_seq.replace('.','')) # Number of none point characters
				if nb_kept == 0:
					# Deletion on full stretch
					if direct_strand:
						# Add information before the stretch
						del_pos = stretch.start - 1
					else:
						# Add information after the stretch
						del_pos = stretch.end + 1
					position_of_deletion(i, del_pos, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, len(del_seq)-nb_kept, allele_seq, pos_list)
	
				elif nb_kept != len(del_seq):
					# Stretch in part deleted
					# Compute positions of substretches (several deletions could be in this stretch)
					del_base_pos = locations_of_substring(del_seq, '.')
					local_stretches = []
					current_substretch = Stretch(False, del_base_pos[0])
					for k in xrange(1,len(del_base_pos)):
						if del_base_pos[k]-1 == del_base_pos[k-1]:
							current_substretch.end += 1
						else:
							current_substretch = Stretch(False, del_base_pos[k])
					# Add the deletion information for each substretch
					for substretch in local_stretches:
						if direct_strand:
							# Add information before the stretch
							del_pos = substretch.start + 1
						else:
							# Add information after the stretch on alignement position, so before in genomics position
							del_pos = substretch.end - 1
						position_of_deletion(i, del_pos, del_pos_list, indel_pos, direct_strand, indel_pos_to_check, align_matrix, len(substretch.pos()), allele_seq, pos_list)

		# Add indels found on deletions positions
		for pos in indel_pos_to_check.keys():
			del_align_matrix[i][del_pos_list.index(pos)] = str(del_align_matrix[i][del_pos_list.index(pos)])+indel_pos_to_check[pos]

	align_matrix = np.c_[align_matrix, del_align_matrix]
	pos_list += del_pos_list

	# Since we added the deletion matrix to the end of the incomplete alignment matrix
	# Positions are not in increasing order anymore, so we will sort it.
	align_matrix = np.transpose(np.transpose(align_matrix)[np.argsort(pos_list)])
	pos_list = np.sort(pos_list)

	# The positions we kept are refering to positions in the alignment and not
	# positions in the reference sequence. The reference sequence can indeed
	# include deletions compared to other sequences. We need to modify the
	# positions list to take this fact into account.

	# List of positions on which a shift occurs.
	# NB : the indices correspond to the number of shifts at this asociated
	# position
	shifts = [0]

	# Find the positions of those shifts
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
			intron_pos_in_file = locations_of_substring(aln_ref_full_seq.replace('.', ''), '|')
			for i in xrange(len(intron_pos_in_file)):
				intron_pos_in_file[i] -= i
			print "Positions of introns in alignment match positions according to ENSEMBL:"
			print intron_pos_in_file == [bio_data[2*k+1] for k in xrange((len(bio_data)-1) / 2)]
			# Check if indels are on the same positions in alignment reference and genome reference
			for i in xrange(len(seq_dict[alignment_ref_seq])):
				if seq_dict[alignment_ref_seq][i] == '.' and seq_dict[genome_ref_seq][i] not in ['.', '*']:
					print "Insertion in genome reference at alignment position " +str(i)+" in file "+file+". Please note that the indels might not be usable in resulting matrix."
				elif seq_dict[genome_ref_seq][i] == '.' and seq_dict[alignment_ref_seq][i] != '.':
					print "Deletion in genome reference at alignment position " +str(i)+" in file "+file+". Please note that the indels might not be usable in resulting matrix."

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