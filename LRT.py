#!/usr/bin/env python

import sys
from Bio.Phylo.PAML.chi2 import cdf_chi2

# This script calculates the p-value for the likelihood ratio test give two
# likelihoods and degrees of freedom.

def usage():
    print "LRT.py <Likelihood 1> <Likelihood 2> <degrees of freedom>"

if len(sys.argv) != 4:
    usage()
    sys.exit()

lnL1 = float(sys.argv[1])
lnL2 = float(sys.argv[2])
degrees = int(sys.argv[3])
LRTstat = 2*(lnL1 - lnL2)
print LRTstat
p = cdf_chi2(degrees, LRTstat)
print p
