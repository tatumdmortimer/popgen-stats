#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict

# This script compares pairwise differences between samples from multisample VCF

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Pairwise Differences')
    parser.add_argument("vcf", help="VCF created from whole genome alignment", 
        type=is_file)
    return parser.parse_args()

def parse_vcf(vcf):
    """ Find variants that are unique to categories """
    variantDict = defaultdict(list)
    inFile = open(vcf, 'r')
    for i,line in enumerate(inFile):
        if i == 2:
            line = line.strip().split()
            strains = line[9:]
        elif i > 2:
            line = line.strip().split()
            variants = line[9:]
            for j,v in enumerate(variants):
                variantDict[strains[j]].append(v)
    inFile.close()
    return variantDict, strains

def pairwise_differences(variantDict, strains):
    with open("pairwiseDifferences.txt", "w") as outfile:
        for i,s in enumerate(strains):
            otherStrains = strains[i+1:]
            for o in otherStrains:
                diff_count = 0
                diff = zip(variantDict[s], variantDict[o])
                for x,y in diff:
                    if x != y:
                        diff_count += 1
                outfile.write("{0}\t{1}\t{2}\n".format(s,o,str(diff_count)))


args = get_args()
variantDict,strains = parse_vcf(args.vcf)
pairwise_differences(variantDict, strains)
