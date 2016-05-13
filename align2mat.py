# -*- coding: utf-8 -*-

################################### Imports ###################################
import numpy as np
import sys


################################## Functions ##################################


##################################### Main ####################################

# Import list of alignment file to process
if len(sys.argv) >= 2 : # List of files file name as an argument
	LIST_FILE = int(sys.argv[1])
else :
	LIST_FILE = "alignments_file_list.txt"

# Put alignment files names in an array (names should not include spaces)
LIST_FILE_DATA = open(LIST_FILE)
FILE_NAMES = LIST_FILE_DATA.readlines()
FILE_NAMES = [name.split() for name in FILE_NAMES]