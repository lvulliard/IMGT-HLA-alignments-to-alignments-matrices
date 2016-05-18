############################################################
	          HLAalignmentToAlignmentMatrix
############################################################
This script aims to convert alignment files, as available
on IMGT-HLA repository 
(ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/)
to alignment matrices (either 0 encoding a similarity to the
reference allele and 1 a difference or 0, 1, 2 and 3
encoding respectively A, T, C and G).
Then, a distance matrix can be computed, by summing the
pairwise differences in the alignment matrix.
############################################################
Run script :
python align2mat.py <source> <codingtype>
python mat2dist.py <source>
With <codingtype> either 0 or 1, for 0 to 3 matrices and 0/1
matrices encoding, and <source> a text file including the 
position of an alignment file on each line. Names should not
include spaces.
############################################################
Output :
.mat - Alignment matrix (binary pickle)
.mat2 - Distance matrix (binary pickle)
.ord - Order of the sequences in the matrices
.pos - Positions kept in matrices
############################################################
CHANGELOG -
v0.1 - Initial comit, the <source> file is read.
v0.2 - Alignment files are read and sequences are cured of
indels and positions with undetermined bases for some
sequences
v0.3 - The alignment matrix is computed, and an additional
script to compute pairwise distances between sequences is
now provided. Matrices are output as pickle files
respectively with a .mat and .mat2 suffixes. The order in
which sequences are reported in the matrices, and the
positions kept to compute the matrices are also presented
in plain text files, respectively with a .ord and .pos
suffixes.
v0.4 - Now the alignment matrix can be encoded as
nucleotides and not only as differences to the reference.