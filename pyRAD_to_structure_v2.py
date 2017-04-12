#!/usr/bin/env python

"""Convert a pyRAD alleles file to a Structure file

author: J. Satler
date: 26 Oct 2016
version: 2

usage: python pyRAD_to_structure_v2.py Infile.alleles
"""

import sys

def generate_loci(alleles):
    """Breaks file into individual loci"""
    with open(alleles, 'r') as a:
        loci_total = []
        locus = 0
        loci = []
        Names = {}
        for line in a:
            line = line.strip()
            if '//' not in line:
                loci.append(line[1:])
                #get names
                sp = line[1:].split()
                if sp[0] not in Names:
                    Names[sp[0]] = []
            else:
                loci_total.append(loci)
                loci = []
                locus += 1
        return Names, loci_total, locus

def haps(locus_lst):
    """convert sequence data to haplotypes,
       treat N's as missing data"""
    poly_loci = []
    mono_loci = []
    for loci in locus_lst:

        l = []
        all = {}
        hap = 1
        haps = []

        #max number of indels at end of sequence
        indel = num_of_indels(loci)
        #print "\nNumber of bps to trim at end: {0}\n\n".format(indel)

        for ind in loci:
            ind = ind.split()

            #trim the right ends off the sequences
            seq_tr = ind[1][:len(ind[1]) - indel]
            #print "O: {0}\nN: {1}\n".format(ind[1], seq_tr)

            #Treat samples with N's as missing data
            if 'N' in seq_tr:
                #print "Sequence contains an N."
                l.append((ind[0], '-9'))
                if '-9' not in haps:
                    haps.append('-9')
                continue

            #first resolved sequence
            if len(all) == 0:
                all[seq_tr] = hap
                l.append((ind[0], hap))
                haps.append(hap)
                hap += 1
            else:
                #haplotype is already recorded
                if all.has_key(seq_tr) == True:
                    l.append((ind[0], all.get(seq_tr)))
                #new haplotype
                else:
                    all[seq_tr] = hap
                    l.append((ind[0], hap))
                    haps.append(hap)
                    hap += 1

        #Remove monomorphic loci
        if len(haps) > 2:
            poly_loci.append(l)
        elif len(haps) == 2 and '-9' not in haps:
            poly_loci.append(l)
        else:
            mono_loci.append(l)
    return poly_loci, mono_loci

def num_of_indels(locus):
    """Count number of indels at end sequence,
        and return number of bps to trim from end"""
    indels = 0
    for ind in locus:
        ind = ind.split()
        #check if sequence ends with gaps
        if ind[1].endswith('-'):
            #Reverse sequence and count number of '-'s
            seq = ind[1][::-1]
            for i in range(len(seq)):
                if seq[i] == '-':
                    continue
                else:
                    #first time we encounter a valid bp
                    if i > indels:
                        indels = i
                        break
                    else:
                        break
    return indels

def add_to_final(samples, haplotypes):
    """Make final list of samples and haplotypes"""
    #keep track of samples with alleles
    sampled = []

    for seq in haplotypes:
        samples[seq[0]].append(seq[1])
        sampled.append(seq[0])
    #if not sequenced, add a '-9' to sample
    for key in samples:
        if key not in sampled:
            samples[key].append("-9")
    return samples

def write_out(str_file):
    """Write out structure file"""
    with open(sys.argv[1] + "_haps.str", 'w') as out:
        #sort samples by name
        for key in sorted(str_file):
            out.write(key + "\t")

            #write out alleles
            alleles = str_file[key]
            for i in range(len(alleles)):
                if i < len(alleles) - 1:
                    out.write(str(alleles[i]) + "\t")
                else:
                    out.write(str(alleles[i]) + "\n")

def Main():
    names, loci, total = generate_loci(sys.argv[1])
    poly, mono = haps(loci)

    print "Polymorphic loci: {0}\tTotal loci: {1}".format(len(poly),
                                                          len(poly)+
                                                          len(mono))
    for locus in poly:
        str = add_to_final(names, locus)
    write_out(str)

if __name__ == "__main__":
    Main()
