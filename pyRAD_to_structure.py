#!/usr/bin/env python

"""Convert a pyRAD alleles file to a Structure file

author: J. Satler
date: 21 Jan 2016

usage: python pyRAD_to_structure.py Infile.alleles
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

        for ind in loci:
            ind = ind.split()
            #Treat samples with N's as missing data
            if 'N' in ind[1]:
                #print "Sequence contains an N."
                l.append((ind[0], '-9'))
                if '-9' not in haps:
                    haps.append('-9')
                continue

            #first resolved sequence
            if len(all) == 0:
                all[ind[1]] = hap
                l.append((ind[0], hap))
                haps.append(hap)
                hap += 1
            else:
                #haplotype is already recorded
                if all.has_key(ind[1]) == True:
                    l.append((ind[0], all.get(ind[1])))
                #new haplotype
                else:
                    all[ind[1]] = hap
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

def haps_with_N(locus_lst):
    """convert sequence data to haplotypes,
       treat N's as missing data"""

    poly_loci = []
    mono_loci = []
    for loci in locus_lst:

        l = []
        all = {}
        hap = 1
        haps = []

        for ind in loci:
            ind = ind.split()
            #Treat samples with N's as missing data
            if 'N' in ind[1]:
                #print "Sequence contains an N."
                l.append((ind[0], '-9'))
                if '-9' not in haps:
                    haps.append('-9')
                continue

            #first resolved sequence
            if len(all) == 0:
                all[ind[1]] = hap
                l.append((ind[0], hap))
                haps.append(hap)
                hap += 1
            else:
                #haplotype is already recorded
                if all.has_key(ind[1]) == True:
                    l.append((ind[0], all.get(ind[1])))
                #new haplotype
                else:
                    all[ind[1]] = hap
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
