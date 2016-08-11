#!/usr/bin/env python

import sys
import getopt

# This script reads in a vcf that has been processed with snpEff.
# Outgroup should be the reference sequence in the vcf.
# The script outputs the synonymous, nonsynonymous, and combined SFS.


def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    vcfFile = None
    numStrains = None
    try:
        opts, args = getopt.getopt(argv, "v:n:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-v':
            vcfFile = arg
        elif opt == '-n':
            numStrains = int(arg)
    return vcfFile, numStrains


def usage():
    print "siteFrequencySpectra.py\n \
        -v <vcf file>\n \
        -n <number of strains (don't include outgroup)>"


def calc_freqs(vcfFile, numStrains):
    nonsynonymous = [0]*numStrains
    synonymous = [0]*numStrains
    intergenic = [0]*numStrains
    non_biallelic = 0
    vcf = open(vcfFile, 'r')
    for line in vcf:
        if line[0] == "#":
            continue
        line = line.strip().split()
        POS = line[1]
        ALT = line[4]
        if len(ALT) > 1:  # skip positions that aren't biallelic
            non_biallelic += 1
            continue
        INFO = line[7]
        outgroup = line[9]
        if outgroup != "0":
            print "VCF in incorrect format."
            print "Outgroup should be reference & first strain in alleles"
            sys.exit()
        alleles = line[10:]
        freq = len(alleles) - alleles.count("0")
        if "synonymous" in INFO:
            synonymous[freq-1] += 1
        elif "missense" in INFO:
            nonsynonymous[freq-1] += 1
        else:
            intergenic[freq-1] += 1
    vcf.close()
    print("{0} SNPs had multiple alternate alleles".format(non_biallelic))
    return synonymous, nonsynonymous, intergenic


def write_outfile(s, ns, ig):
    outfile = open("sfs.txt", "w")
    outfile.write("Frequency\tSynonymous\tNonsynonymous\tIntergenic\tCombined\n")
    for i in range(len(s)):
        outfile.write("%i\t%i\t%i\t%i\t%i\n" % (i+1,
                                                s[i],
                                                ns[i],
                                                ig[i],
                                                s[i] + ns[i] + ig[i]))


vcfFile, numStrains = get_arguments(sys.argv[1:])

if vcfFile is None or numStrains is None:
    usage()
    sys.exit()

synonymous, nonsynonymous,intergenic = calc_freqs(vcfFile, numStrains)
write_outfile(synonymous, nonsynonymous, intergenic)
