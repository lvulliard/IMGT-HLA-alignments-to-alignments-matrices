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
Basic usage:
align2mat.py --help
mat2dist.py --help
align2mat.py --list <source> --coord <grch>
mat2dist.py --list <source>
With <grch> the version of the genome in which you want the
coordinates of your alignment matrices, and <source> a text
file including the path to an alignment file on each line.
Names should not include spaces.
Option --help show all the other options available.
############################################################
Output:
.mat - Alignment matrix (binary pickle)
.mat2 - Distance matrix (binary pickle)
.ord - Order of the sequences in the matrices
.pos - Positions kept in matrices
############################################################
Additional information:
In order to convert positions from alignment coordinates to
genomic coordinates, you need a file describing for each 
alignment file the structure of present features.
For instance, if you want to process the alignment file
"X.txt" and get the positions as found in GRCh38, you will 
need a file called "X.txt.GRCh38.bio". On the first line of
this file, you will have either 0 or 1, respectively if the
corresponding gene is coding on the reverse strand or on the
direct strand. The second line corresponds to the start of
the first base of the reference sequence in the alignment.
Then odd lines correspond to number of bases until the next
intron, and even lines to the length of an intron. The last
line corresponds to the position of the last base of the
reference sequence. All that information can be found on 
ENSEMBL database.
Finally, to use this information you will need to add -c 38
or --coord 38 when calling align2mat.py.
You can check if positions of introns in the alignments,
encoded by pipe characters, are matching expectations
according to the .GRChX.bio file by adding the -a or the
--addchecks tag.
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
v0.5 - Positions are now corrected to take into account
deletions in the reference genome
v0.6 - Improved how user inputs script options