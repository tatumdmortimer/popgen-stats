#!/usr/bin/env python

import sys
import os
import argparse
import numpy 
from numpy import array

# This script randomly sames a site frequency spectrum with replacement and 
# creates input files for dadi and prfreq

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Sample SFS with replacement')
    parser.add_argument("sfs", help="SFS file to be sampled",
                    type=argparse.FileType('r'))
    parser.add_argument("n", help="Number of samples", type=int)
    parser.add_argument("-p", "--prfreq", help="Create prfreq input files",
                        action="store_true")
    parser.add_argument("-d", "--dadi", help="Create dadi input files",
                        action="store_true")
    return parser.parse_args()
    
args = get_args()

freqs = []
num = 0
# read in SFS
for i,line in enumerate(args.sfs):
    line = line.strip()
    freqs += [i+1]*int(line)
    num = i + 1

# create directories for output files
if args.dadi:
    os.mkdir("dadi_samples")
if args.prfreq:
    os.mkdir("prfreq_samples")
    
# randomly sample SFS
freqs = array(freqs)
for i in range(args.n):
    sample = numpy.random.choice(freqs, len(freqs))
    sample_sfs = numpy.bincount(sample)
    if args.dadi:
        dadi_out = open("dadi_samples/sample%i.txt" % (i+1), "w")
        dadi_out.write("# SFS sample %i\n" % (i+1))
        dadi_out.write("# Outgroup: Troy\n")
        dadi_out.write("%i unfolded\n" % (num + 1))
        dadi_out.write(" ".join(str(k) for k in sample_sfs) + "\n")
        dadi_out.close()
    if args.prfreq:
        prfreq_out = open("prfreq_samples/sample%i.txt" % (i+1), "w")
        for m in sample_sfs[1:-1]:
            prfreq_out.write(str(m) + "\n")
        prfreq_out.close()
