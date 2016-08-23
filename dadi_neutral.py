#!/usr/bin/env python

import sys
import argparse
import numpy
from numpy import array
import dadi
from datetime import datetime


def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Run dadi analysis')
    parser.add_argument("sfs", help="SFS in dadi format",
                    type=argparse.FileType('r'))
    parser.add_argument("-f", "--fold", help="Fold the SFS before analysis",
                        action="store_true")
    return parser.parse_args()

args = get_args()

startTime = datetime.now()

data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()

# Neutral model
pts_l = [110,120,130] # grid point settings
neutral_func = dadi.Demographics1D.snm
neutral_params = array([])
neutral_upper_bound = []
neutral_func_ex = dadi.Numerics.make_extrap_log_func(neutral_func)
neutral_model = neutral_func_ex(neutral_params, ns, pts_l)

neutral_ll = dadi.Inference.ll_multinom(neutral_model, data)

print "Neutral model log-likelihood: %f" % neutral_ll

# Output SFS for data
data_sfs_file = open("observedSFS.txt", "w")
for i in range(1,len(data)-1):
    data_sfs_file.write(str(data[i]) + '\n')
data_sfs_file.close()

# Output SFS for neutral model
neutral_sfs = dadi.Inference.optimally_scaled_sfs(neutral_model, data)
neutral_sfs_file = open("neutralModelSFS.txt", 'w')
for i in range(1,len(neutral_sfs)-1):
    neutral_sfs_file.write(str(neutral_sfs[i]) + '\n')
neutral_sfs_file.close()

print datetime.now() - startTime
