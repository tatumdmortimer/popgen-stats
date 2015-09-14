#!/usr/bin/env python

import sys
import os
import argparse
import glob
import egglib

# This script reads in a fasta alignment or a directory of alignments 
# and a list of outgroup sequences. Alignments in directory must end in .fasta.
# Alignments must also be in frame coding sequence. 
# It outputs population genetic statistics for including theta, pi,
# piN, piS, Tajima's D, and MK test table for each outgroup provided.
# This script requires egglib installed with the Bio++ libraries


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

def get_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Calculate diversity and selection statistic")
    input_files = parser.add_mutually_exclusive_group(required=True)
    input_files.add_argument('-a', '--alignment', help = 'Alignment to calculate statistics',
        type=is_file)
    input_files.add_argument('-d', '--directory', help = 'Directory of alignments',
        type=is_dir)
    parser.add_argument('-f', '--frame', help = 'Alignment is in correct reading frame', action='store_true')
    parser.add_argument('-o', '--outgroup', help = 'Outgroup(s) for McDonald-Kreitman Test', 
        type=str, nargs='+')
    return parser.parse_args()
     
def calc_stats(alignment, outgroup):
    statDict = {}
    a = egglib.Align(alignment)
    for i in range(a.ns()):
        a.sequence(i, sequence=a.sequence(i).upper())
    polyDict = a.polymorphism()
    statDict['theta'] = polyDict['thetaW']
    statDict['pi'] = polyDict['Pi']
    statDict['tajimaD'] = polyDict['D']
    if args.frame:
        if len(a.sequence(1))%3 != 0:
            print("The following alignment is not in frame:")
            print(alignment)
            sys.exit()
        polyDictBPP = a.polymorphismBPP(dataType=4)
        statDict['piN'] = polyDictBPP['PiNS']
        statDict['piS'] = polyDictBPP['PiS']
    if outgroup is not None:
        for o in outgroup:
            temp_a = a.extract(0, len(a.sequence(1)) - 3) # remove stop codon
            otherOutgroups = outgroup[:]
            otherOutgroups.remove(o)
            for otherOutgroup in otherOutgroups: 
                index = temp_a.find(otherOutgroup, strict=False)
                if index != None:
                    del temp_a[index]
            try:
                temp_a.group(temp_a.find(o, strict = False), group=999)
            except IndexError:
                print("The following outgroup is not present in alignment")
                print(o, a.find(o, strict = False))
                sys.exit()
            polyDictBPP = temp_a.polymorphismBPP(dataType=4)
            statDict['MK_'+o] = polyDictBPP['MK']
            statDict['NI_'+o] = polyDictBPP['NI']
    return statDict

def write_outfile(alignDict, outgroup):
    outfile = open('selectionStats.txt', 'w')
    outfile.write('Alignment\tTheta\tPi\tTajimasD')
    if args.frame:
        outfile.write('\tPiN\tPiS')
        if outgroup is not None:
            for o in outgroup:
                outfile.write('\tMK_' + o)
                outfile.write('\tNI_' + o)
    outfile.write('\n')
    for a in alignDict:
        s = alignDict[a]
        if len(s) == 0:
            print a, " is not in frame"
            continue
        outfile.write('%s\t%s\t%s\t%s' % (a, s['theta'],s['pi'], s['tajimaD']))
        if args.frame:
            outfile.write('\t%s\t%s' % (s['piN'], s['piS']))
            if outgroup is not None:
                for o in outgroup:
                    outfile.write('\t' + str(s['MK_' + o]))
                    outfile.write('\t' + str(s['NI_' + o]))
        outfile.write('\n')
    outfile.close()

args = get_arguments()

alignDict = {}
# Check if alignment or directory was given and calculate stats accordingly
if args.alignment is None:
    for align in glob.glob(args.directory + '*.fasta'):
        alignName = os.path.splitext(align)[0].replace(directory, "")
        alignDict[alignName] = calc_stats(align, args.outgroup)

else:
    alignName = os.path.splitext(args.alignment)[0]
    alignDict[alignName] = calc_stats(args.alignment, args.outgroup)

write_outfile(alignDict, args.outgroup)
