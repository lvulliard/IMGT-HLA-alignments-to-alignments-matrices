############################################################
	          HLAalignmentToAlignmentMatrix
############################################################
This script aims to convert alignment files, as available
on IMGT-HLA repository 
(ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/)
to alignment matrices (0 encoding a similarity to the
reference allele, and 1 a difference).
############################################################
Run script : python align2mat.py <source>
With <source> a text file including the position of an
alignment file on each line. Names should not include
spaces.
############################################################
CHANGELOG -
v0.1 - Initial comit, the <source> file is read.
v0.2 - Alignment files are read and sequences are cured of
indels and positions with undetermined bases for some
sequences