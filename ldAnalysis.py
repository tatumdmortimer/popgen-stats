#!/usr/bin/env python

import sys
import os
import argparse
import egglib
import subprocess

# This script uses the egglib package to calculate LD statistics within and
# between genes in a directory genes

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

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
    parser.add_argument("directory", help="Directory of alignments", 
        type=is_dir, action=FullPaths)
    return parser.parse_args()

def interGeneLD(directory):
    files = listdir_fullpath(directory)
    if not os.path.exists("interLD"):
        os.makedirs("interLD")
    for f in range(len(files)):
        f_name = files[f]
        f_gene = os.path.basename(f_name).split("_")[0]
        for g in range(f + 1, len(files)):
            g_name = files[g]
            g_gene = os.path.basename(g_name).split("_")[0]
            subprocess.call(["/opt/PepPrograms/egglib-py-2.1.7/scripts/egglib",
                "interLD",
                "align1={0}".format(f_name),
                "align2={0}".format(g_name),
                "output=interLD/{0}_{1}_LDout.txt".format(f_gene, g_gene)])
        

def withinGeneLD(directory):
    outfile = open("LDstats.txt", "w")
    outfile.write("Gene\tPosition1\tPosition2\tDistance\tD'\tR2\n")
    for f in listdir_fullpath(directory):
        a = egglib.Align(f)
        gene = os.path.basename(f).split("_")[0]
        siteIndices = a.polymorphism()['siteIndices']
        ldStats = a.matrixLD()
        for i in range(len(siteIndices)):
            for j in range(i+1, len(siteIndices)):
                pos1 = siteIndices[i]
                pos2 = siteIndices[j]
                distance = ldStats['d'][pos1][pos2]
                dPrime = ldStats['Dp'][pos1][pos2]
                rSquared = ldStats['r2'][pos1][pos2]
                outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    gene, pos1, pos2, distance, dPrime, rSquared))
    outfile.close()


args = get_args()
withinGeneLD(args.directory)
interGeneLD(args.directory)
