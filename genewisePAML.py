#!/usr/bin/env python

import sys
import os
import getopt
import glob
import egglib
import egglib.wrappers as wrappers
from Bio.Phylo.PAML.chi2 import cdf_chi2

# This script reads in a fasta alignment or a directory of alignments.
# Alignments in directory must end in .fasta.
# Alignments must also be in frame coding sequence. 

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignmentFile = None
    alignmentDirectory = None
    try:
        opts, args = getopt.getopt(argv, "a:d:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignmentFile = arg
        elif opt == '-d':
            alignmentDirectory = arg
    return (alignmentFile, alignmentDirectory)

def usage():
    print "genewisePAML.py\n \
        -a <fasta alignment>\n \
        -d <directory of fasta alignments>"

def calc_tree(alignment):
    a = egglib.Align(alignment)
    a = a.extract(0, len(a.sequence(1)) - 3) # remove stop codon
    if len(a.sequence(1))%3 != 0:
        print(alignment, " not in frame")
        return (a, None)
    if a.ns() < 3:
        return (a, None)
    for i in range(a.ns()):
        a.sequence(i, sequence=a.sequence(i).upper())
    tree, loglk = wrappers.phyml(a)
    return (a, tree)

def run_paml(a, tree, alignName, outfile, neutralFile):
    try:
        codemlInstance = wrappers.Codeml(a, tree)
        neutral = codemlInstance.fit("M1a")
        positive = codemlInstance.fit("M2a")
    except ValueError, e:
        print (alignName, e)
        return
    fitDict = {}
    n = True
    LRTstat = 2*(positive["lnL"] - neutral["lnL"])
    if LRTstat > 0:
        p = cdf_chi2(2, LRTstat)
        if p < 0.05:
            n = False
        if not n:
            outfile.write("%s\t%f\n" % (alignName, p))
        else:
            neutralFile.write("%s\t%f\n" % (alignName, p))
       
    
alignment, directory = get_arguments(sys.argv[1:])
alignDict = {}
# Check if alignment or directory was given and calculate stats accordingly
outfile = open("pamlResults_significant.txt", "w")
neutralFile = open("pamlResults_neutral.txt", "w")
if alignment is None:
    if directory is None:
        usage()
        sys.exit()
    else:
        if not os.path.isdir(directory):
            print("This is not a valid directory")
            usage()
            sys.exit()
        totalAlign = len(glob.glob(directory + '*.fasta'))
        for i,align in enumerate(glob.glob(directory + '*.fasta')):
            alignName = os.path.splitext(align)[0].replace(directory, "")
            a, tree = calc_tree(align)
            if tree is None:
                continue
            run_paml(a, tree, alignName, outfile, neutralFile)
            percentageDone = 100*(float(i)/totalAlign)
            if i%10 >= 0 and i%10 <= 1:
                print "Analysis is %f percent complete" % (percentageDone)
                outfile.flush()
                os.fsync(outfile)
                neutralFile.flush()
                os.fsync(neutralFile)

elif alignment is not None:
    if directory is not None:
        print "Must only input an alignment or a directory"
        usage()
        sys.exit()
    else:
        if not os.path.isfile(alignment):
            print "The file does not exist."
            usage()
            sys.exit()
        alignName = os.path.splitext(alignment)[0]
        a, tree = calc_tree(alignment)
        if tree is None:
            print ("This alignment does not have a tree")
        run_paml(a, tree, alignName, outfile, neutralFile)

outfile.close()
neutralFile.close()
