# -*- coding: utf-8 -*-

################################### Imports ###################################
import numpy as np
import sys, pickle


################################## Functions ##################################

##################################### Main ####################################

# Import list of alignment file to process
if len(sys.argv) >= 2 : # List of files file name as an argument
	LIST_FILE = sys.argv[1].split()[0]
else :
	LIST_FILE = "alignments_file_list.txt"

# Put alignment files names in an array (names should not include spaces)
LIST_FILE_DATA = open(LIST_FILE)
FILE_NAMES = (name.split()[0] for name in LIST_FILE_DATA.readlines())

# Loop on every file to process
for file in FILE_NAMES:
	print "Processing "+file

	# Read associated matrix
	align_matrix = pickle.load(open(file+".mat",'rb'))
	# Number of sequences
	nbseq = len(align_matrix)
	# Number of kept positions
	nbpos = len(align_matrix[0])
	# Create distance matrix
	dist_matrix = np.empty([nbseq, nbseq])
	# Compute pairwise number of differences in the alignment matrix
	for i in xrange(nbseq):
		j = 0
		while j < i:
			dist_ij = 0
			for k in xrange(nbpos):
				if align_matrix[i][k] != align_matrix[j][k]:
					dist_ij += 1
			dist_matrix[i][j] = dist_ij
			dist_matrix[j][i] = dist_ij
			j += 1
		dist_matrix[i][i] = 0

	# Store distances in pickle .mat2 file
	pickle.dump(dist_matrix, open(file+".mat2",'wb'))