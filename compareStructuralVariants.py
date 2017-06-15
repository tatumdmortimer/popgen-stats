#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict

# This script filters VCFs output by pindel and compares structural varaiants
# between genomes of groups input by the user

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
    parser.add_argument("categories", help="File describing genome categories",
        type=is_file)
    return parser.parse_args()

def read_VCFs(strain):
    variantDict = {}
    file_endings = ["_D.recode.vcf", "_INV.recode.vcf", "_SI.recode.vcf",
    "_TD.recode.vcf"]
    for ext in file_endings:
        with open(strain + ext, "r") as infile:
            for line in infile:
                if line[0] != "#":
                    line = line.strip().split()
                    position = line[1]
                    alt = line[4]
                    variantDict[position] = alt
    return variantDict


def read_cat_file(genomeCatFile):
    """ Read in genome categories and create dictionary of category name and 
    genomes in that category. Additionally call function to read VCFs."""
    inFile = open(genomeCatFile, 'r')
    catDict = defaultdict(list)
    strain_variant_dict = {}
    for line in inFile:
        line = line.strip()
        entries = line.split()
        genome = entries[0]
        cat = entries[1]
        catDict[cat].append(genome)
        strain_variant_dict[genome] = read_VCFs(genome)
    inFile.close()
    return (catDict, strain_variant_dict)

def analyze_variants(catDict, strain_variant_dict):
    """ Find variants that are unique to categories """
    sharedVariantDict = {}
    allVariantDict = {}
    ##Find variants that are common to all strains in category
    for c in catDict:
        for i,g in enumerate(catDict[c]):
            var = set([(k,v) for k, v in strain_variant_dict[g].iteritems()])
            if i == 0:
                sharedVariantDict[c] = var
                allVariantDict[c] = var
            else:
                sharedVariantDict[c] = sharedVariantDict[c] & var
                allVariantDict[c] = allVariantDict[c] | var
    ##Find variants that are unique to each category
    uniqueVariantDict = {}
    for c in catDict:
        unique = sharedVariantDict[c].copy()
        for otherCat in catDict:
            if otherCat != c:
                unique -= allVariantDict[otherCat]
        uniqueVariantDict[c] = unique
    return sharedVariantDict, uniqueVariantDict

def write_variants(sharedVariantDict, uniqueVariantDict):
    for cat in sharedVariantDict:
        outFile = open(cat + "_sharedVariants.txt", "w")
        for pos,alt in sharedVariantDict[cat]:
            outFile.write("{0}\t{1}\n".format(pos,alt))
        outFile.close()
        uniqueFile = open(cat + "_uniqueVariants.txt", "w")
        for pos,alt in uniqueVariantDict[cat]:
            uniqueFile.write("{0}\t{1}\n".format(pos,alt))
        uniqueFile.close()

args = get_args()

catDict, strainVariantDict = read_cat_file(args.categories)
sharedVariantDict, uniqueVariantDict = analyze_variants(catDict, 
    strainVariantDict)
write_variants(sharedVariantDict, uniqueVariantDict)
