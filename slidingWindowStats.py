#!/usr/bin/env python

import sys
import os
import getopt
import egglib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# This script reads in an alignment and calculates diversity and selection
# statistics based on the window width and window step given by the user.
# Will calculate Fay and Wu's H if an outgroup is provided.

def get_arguments(argv):
    if len(argv) == 0:
        usage()
        sys.exit(2)
    alignment = None
    winWidth = 1000
    winStep = 300
    outgroup = None
    try:
        opts, args = getopt.getopt(argv, "a:w:s:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-a':
            alignment = arg
        elif opt == '-w':
            winWidth = int(arg) 
        elif opt == '-s':
            winStep = int(arg)
        elif opt == '-o':
            outgroup = arg
    return (alignment, winWidth, winStep, outgroup)

def usage():
    print "slidingWindowStats.py\n \
        -a <fasta alignment>\n \
        -w <window width default = 1000>\n \
        -s <window step default = 300>\n \
        -o <outgroup>"

def calc_stats(a):
    statDict = {}
    polyDict = a.polymorphism()
    statDict['theta'] = polyDict['thetaW']
    statDict['pi'] = polyDict['Pi']
    statDict['tajimaD'] = polyDict['D']
    statDict['FayWuH'] = polyDict['H']
    return statDict

alignment, winWidth, winStep, outgroup = get_arguments(sys.argv[1:])

if alignment is None:
    usage()
    sys.exit()

outfile = open('windowStats_' + os.path.splitext(alignment)[0] + '.txt', 'w')
outfile.write("Start\tStop\tTheta\tPi\tTajimasD\tFay&WuH\n")
align = egglib.Align(alignment)
for i in range(align.ns()):
    align.sequence(i, sequence=align.sequence(i).upper())
if outgroup is not None:
    align.group(align.find(outgroup, strict = False), group = 999)
start = 0
stop = winWidth

location = []
TD = []

for window in align.slider(winWidth, winStep):
    stats = calc_stats(window)
    start += winStep
    stop += winStep
    outfile.write("%i\t%i\t%s\t%s\t%s\t%s\n" % (start, stop, stats['theta'],
    stats['pi'], stats['tajimaD'], stats['FayWuH']))
    location.append((start + stop)/2)
    TD.append(stats['tajimaD'])
outfile.close()

plt.plot(location, TD)
plt.xlabel('Location')
plt.ylabel('Tajima\'s D')
plt.savefig("slidingWindowTajimasD.png")
plt.close()
