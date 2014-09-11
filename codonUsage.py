#!/usr/bin/env python

import sys
import os
import glob
import subprocess
from Bio import SeqIO

# This script computes CAI (codon adaptation index) for a directory of genes
# Files in directory must end in .fasta.
# Must also be in frame coding sequence, but doesn't need to be aligned..
# It outputs population genetic statistics for including theta, pi,
# A file containing the preferred codon usage table created by emboss cusp
# should also be provided.                                      

def usage():
    print "codonUsage.py <directory of fasta alignments> <codon usage table>"


# Check if directory and codon usage table were given
if len(sys.argv) != 3:
    usage()
    sys.exit()



directory = sys.argv[1]
cut = sys.argv[2]

if not os.path.isdir(directory):
    print "Directory of alignments not found"
    usage()
    sys.exit()

if not os.path.exists(cut):
    print "Codon usage table file not found"
    usage()
    sys.exit()

geneDict = {}
strainList = []
files = glob.glob(directory + "*.fasta")
os.mkdir("tmp") # make temporary directory to store intermediate files
for f in files:
    # calculate CAI for each sequence in the file using EMBOSS program cai
    outfile = "tmp/" + f.split("/")[-1].split(".")[0] + ".cai"
    subprocess.call(["/home/peplab/src/EMBOSS-6.5.7/emboss/cai",
                        "-seqall", f, "-cfile", cut,
                        "-outfile", outfile]) 
caiFiles = glob.glob("tmp/*.cai")

for i, caiFile in enumerate(caiFiles):
    caiDict = {}
    gene = caiFile.split("/")[1].split(".")[0]
    cf = open(caiFile, "r")
    for line in cf:
        line = line.strip().split()
        strain = line[1].split("_")[0]
        if i == 0:
            strainList.append(strain)
        CAI = line[3]
        caiDict[strain] = CAI
    geneDict[gene] = caiDict
    cf.close()
    os.remove(caiFile)

outFile = open("CAI.txt", "w")

outFile.write("Gene\t" + "\t".join(strainList) + "\n")
genes = geneDict.keys()
for g in genes:
    outFile.write(g)
    for s in strainList:
        outFile.write("\t" + geneDict[g][s])
    outFile.write("\n")

outFile.close()

os.rmdir("tmp/")
