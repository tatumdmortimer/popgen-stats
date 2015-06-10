#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict

# This script filters protein groups output from OrthoMCL and compares the core
# genomes of groups input by the user

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Compare core genomes')
    parser.add_argument("vcf", help="VCF created from whole genome alignment", 
        type=is_file)
    parser.add_argument("categories", help="File describing genome categories",
        type=is_file)
    return parser.parse_args()

def read_cat_file(genomeCatFile):
    """ Read in genome categories and create dictionary of category name and 
    genomes in that category"""
    inFile = open(genomeCatFile, 'r')
    catDict = {}
    catSet = set()
    for line in inFile:
        line = line.strip()
        entries = line.split()
        genome = entries[0]
        cat = entries[1]
        catDict[genome] = cat
        catSet.add(cat)
    inFile.close()
    return (catDict, catSet)

def analyze_variants(catDict, catSet, vcf):
    """ Find variants that are unique to categories """
    variantDict = defaultdict(list)
    inFile = open(vcf, 'r')
    for i,line in enumerate(inFile):
        if i == 2:
            line = line.strip().split()
            strains = line[9:]
            categories = [catDict[strain] for strain in strains]
        elif i > 2:
            tempDict = defaultdict(set)
            line = line.strip().split()
            position = line[1]
            variants = line[9:]
            for j,v in enumerate(variants):
                tempDict[categories[j]].add(v)
            for c in tempDict:
                if len(tempDict[c]) != 1:
                    continue
                for otherCat in tempDict:
                    if otherCat == c:
                        continue
                    elif tempDict[c].issubset(tempDict[otherCat]):
                        break
                    
                else:
                    variantDict[c].append(position)

    inFile.close()
    return variantDict

def write_variants(variantDict):
    for cat in variantDict:
        outFile = open(cat + "_variants.txt", "w")
        for variant in variantDict[cat]:
            outFile.write(variant + "\n")
        outFile.close()

args = get_args()

catDict,catSet = read_cat_file(args.categories)
variantDict = analyze_variants(catDict, catSet, args.vcf)
write_variants(variantDict)
