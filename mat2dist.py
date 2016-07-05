#!/usr/bin/env python
# -*- coding: utf-8 -*-
##################################### Help ####################################
"""Convert IMGT-HLA alignments to alignment matrices.

Usage:
  align2mat.py [-l <path>]
  align2mat.py --help
  align2mat.py --version

Options:
  -l <path>, --list <path>    Path to a list of files to process. [default: alignments_files_list.txt]         
  -h --help                   Show this screen.
  -v --version                Show version.

"""

################################### Imports ###################################
from docopt import docopt
import numpy as np
import sys, pickle


################################## Functions ##################################

##################################### Main ####################################

# Import user input
arguments = docopt(__doc__, version='v0.6')
LIST_FILE = arguments["--list"] # List of files file name

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