#!/usr/bin/env python

import sys
import os
import argparse
import random

# This script reads an output file from slim and creates input for dadi 

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Converts slim output to dadi\
 input')
    parser.add_argument("slim", help="slim output file",
                    type=argparse.FileType('r'))
    parser.add_argument("dadi", help="name for dadi input file",
                    type=argparse.FileType('w'))
    return parser.parse_args()
    
def create_sfs_slim(slimFile):
    pop_section = False
    mut_section = False
    sfs = None
    for line in slimFile:
        if line.strip() == "Populations:":
            pop_section = True
            continue
        elif line.strip() == "Mutations:":
            pop_section = False
            mut_section = True
            continue
        elif line.strip() == "Genomes:":
            mut_section = False
            break
        elif pop_section:
            line = line.strip().split()
            n = int(line[1])*2
            sfs = [0]*(n+1)
        elif mut_section:
            line = line.strip().split()
            freq = int(line[7])
            sfs[freq] += 1
    return sfs

def write_dadi(sfs, dadi_out):
    dadi_out.write("# SFS from slim simulation\n")
    dadi_out.write("%i unfolded\n" % (len(sfs)))
    dadi_out.write(" ".join(str(k) for k in sfs) + "\n")
    dadi_out.close()

args = get_args()

sfs = create_sfs_slim(args.slim)
write_dadi(sfs, args.dadi)
