#!/usr/bin/env python

import sys
import os
import getopt
import glob
import egglib

# This script reads in a fasta alignment or a directory of alignments 
# and a list of outgroup sequences. Alignments in directory must end in .fasta.
# Alignments must also be in frame coding sequence. 
# It outputs population genetic statistics for including theta, pi,
# piN, piS, Tajima's D, and MK test table for each outgroup provided.
# This script requires egglib installed with the Bio++ libraries

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignmentFile = None
    alignmentDirectory = None
    outgroups = None
    try:
        opts, args = getopt.getopt(argv, "a:d:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignmentFile = arg
        elif opt == '-d':
            alignmentDirectory = arg
        elif opt == '-o':
            outgroups = arg
    return (alignmentFile, alignmentDirectory, outgroups)

def usage():
    print "FilterOrthoMCLGroups.py\n \
        -a <fasta alignment>\n \
        -d <directory of fasta alignments>\n \
        -o <'list, of, outgroups'>"

def calc_stats(alignment, og, outgroup):
    print alignment
    statDict = {}
    if og:
        for o in outgroup:
            a = egglib.Align(alignment)
            otherOutgroups = outgroup[:]
            otherOutgroups.remove(o)
            for otherOutgroup in otherOutgroups: 
                index = a.find(otherOutgroup, strict=False)
                if index != None:
                    del a[index]
            a.group(a.find(o, strict = False), group=999)
            for i in range(a.ns()):
                a.sequence(i, sequence=a.sequence(i).upper())
            polyDict = a.polymorphism()
            polyDictBPP = a.polymorphismBPP(dataType=4)
            statDict['theta'] = polyDict['thetaW']
            statDict['pi'] = polyDict['Pi']
            statDict['tajimaD'] = polyDict['D']
            statDict['piN'] = polyDictBPP['PiNS']
            statDict['piS'] = polyDictBPP['PiS']
            statDict['MK_'+o] = polyDictBPP['MK']
    else:
        a = egglib.Align(alignment)
        for i in range(a.ns()):
             a.sequence(i, sequence=a.sequence(i).upper())
        polyDict = a.polymorphism()
        polyDictBPP = a.polymorphismBPP(dataType=4)
        statDict['theta'] = polyDict['thetaW']
        statDict['pi'] = polyDict['Pi']
        statDict['tajimaD'] = polyDict['D']
        statDict['piN'] = polyDictBPP['PiNS']
        statDict['piS'] = polyDictBPP['PiS']
       
    return statDict

def write_outfile(alignDict, og, outgroup):
    outfile = open('selectionStats.txt', 'w')
    outfile.write('Alignment\tTheta\tPi\tPiN\tPiS\tTajimasD')
    if og:
        for o in outgroup:
            outfile.write('\tMK_' + o)
    outfile.write('\n')
    for a in alignDict:
        s = alignDict[a]
        outfile.write('%s\t%s\t%s\t%s\t%s\t%s' % (a, s['theta'],s['pi'],s['piN'],
                                                s['piS'],s['tajimaD']))
        if og:
            for o in outgroup:
                outfile.write('\t' + str(s['MK_' + o]))
        outfile.write('\n')
    outfile.close()

alignment, directory, outgroup = get_arguments(sys.argv[1:])

print outgroup
# Check if there are outgroups
og = False
if outgroup is not None:
    og = True
    outgroup = outgroup.split(', ')

alignDict = {}
i = 0
# Check if alignment or directory was given and calculate stats accordingly
if alignment is None:
    if directory is None:
        usage()
        sys.exit()
    else:
        for align in glob.glob(directory + '*.fasta'):
            print i
            i += 1
            alignName = os.path.splitext(align)[0].replace(directory, "")
            alignDict[alignName] = calc_stats(align, og, outgroup)

elif alignment is not None:
    if directory is not None:
        print "Must only input an alignment or a directory"
        usage()
        sys.exit()
    else:
        alignName = os.path.splitext(alignment)[0]
        alignDict[alignName] = calc_stats(alignment, og, outgroup)

write_outfile(alignDict, og, outgroup)
