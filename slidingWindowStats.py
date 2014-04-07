#!/usr/bin/env python

import sys
import os
import getopt
import egglib

# This script reads in an alignment and calculates diversity and selection
# statistics based on the window width and window step given by the user.

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignment = None
    winWidth = 1000
    winStep = 300
    try:
        opts, args = getopt.getopt(argv, "a:w:s:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignment = arg
        elif opt == '-w':
            winWidth = int(arg) 
        elif opt == '-o':
            winStep = int(arg)
    return (alignment, winWidth, winStep)

def usage():
    print "slidingWindowStats.py\n \
        -a <fasta alignment>\n \
        -w <window width default = 1000>\n \
        -s <window step default = 300>"

def calc_stats(a):
    statDict = {}
    polyDict = a.polymorphism()
    statDict['theta'] = polyDict['thetaW']
    statDict['pi'] = polyDict['Pi']
    statDict['tajimaD'] = polyDict['D']
    return statDict

alignment, winWidth, winStep = get_arguments(sys.argv[1:])

if alignment is None:
    usage()
    sys.exit()

outfile = open('windowStats_' + os.path.splitext(alignment)[0] + '.txt', 'w')
outfile.write("Start\tStop\tTheta\tPi\tTajimasD\n")
align = egglib.Align(alignment)
for i in range(align.ns()):
    align.sequence(i, sequence=align.sequence(i).upper())
start = 0
stop = winWidth
for window in align.slider(winWidth, winStep):
    stats = calc_stats(window)
    start += winStep
    stop += winStep
    outfile.write("%i\t%i\t%s\t%s\t%s\n" % (start, stop, stats['theta'],
    stats['pi'], stats['tajimaD']))
outfile.close()
