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
	# Read file, skipping header
	file_data = open(file)
	lines = file_data.readlines()[6:]

	# Parse file
	line_index = 0
	while(line_index < len(lines)):
		# Go to next line
		line_index +=1