#pyRAD_to_Structure_v2.py

This script will take pyRAD output (*.alleles file) and convert it into an
input file for STRUCTURE. It uses full haplotypes, and if an 'N' is present
in the sequence, the haplotype is treated as missing ('-9'). This is conservative
but should only impact a small amount of data. In addition, missing data at the end
of the sequences ('-') are not treated as characters, so are ignored when determining
the total number of haplotypes (at a locus) and how they are assigned to individuals.