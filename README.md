#pyRAD_to_Structure

This script will take pyRAD output (*.alleles file) and convert it into an
input file for STRUCTURE. It uses the full haplotypes, and if an 'N' is present
in the sequence, the haplotype is treated as missing ('-9'). This is conservative
but should only impact a small amount of data.