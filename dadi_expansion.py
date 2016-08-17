#!/usr/bin/env python

import sys
import argparse
from Bio.Phylo.PAML.chi2 import cdf_chi2
import numpy
from numpy import array
import dadi
from math import *
import cmath
from datetime import datetime


def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Run dadi analysis')
    parser.add_argument("sfs", help="SFS in dadi format",
                        type=argparse.FileType('r'))
    parser.add_argument("-f", "--fold", help="Fold the SFS before analysis",
                        action="store_true")
    return parser.parse_args()


def rep_ll_expansion(data, ns, expansion_params, expansion_upper_bound,
                     expansion_lower_bound, expansion_func_ex, pts_l):
    expansion_vals = []
    outfile = open("dadi_expansionModel.txt", "w")
    outfile.write("nu\ttau\tll\n")
    for i in range(100):
        expansion_p0 = dadi.Misc.perturb_params(
            expansion_params, fold=1, upper_bound=expansion_upper_bound)

        expansion_popt = dadi.Inference.optimize_log(
            expansion_p0, data, expansion_func_ex, pts_l,
            lower_bound=expansion_lower_bound,
            upper_bound=expansion_upper_bound, maxiter=100)
        expansion_model = expansion_func_ex(expansion_popt, ns, pts_l)
        expansion_ll_opt = dadi.Inference.ll_multinom(expansion_model, data)
        outfile.write(
            "{0}\t{1}\t{2}\n".format(
                expansion_popt[0],
                expansion_popt[1],
                expansion_ll_opt))
        expansion_vals.append(expansion_ll_opt)
    avg_ll_expansion = sum(expansion_vals) / float(len(expansion_vals))
    outfile.close()
    return avg_ll_expansion


startTime = datetime.now()

pts_l = [110, 120, 130]
expansion_func = dadi.Demographics1D.two_epoch
expansion_params = array([2, 0.05])
expansion_upper_bound = [100, 10]
expansion_lower_bound = [1e-2, 0]
expansion_func_ex = dadi.Numerics.make_extrap_log_func(expansion_func)
k = len(expansion_params)

args = get_args()
data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()

avg_ll_expansion = rep_ll_expansion(data, ns, expansion_params, expansion_upper_bound,
                                    expansion_lower_bound, expansion_func_ex, pts_l)
print "Average optimized log-likelihood:", avg_ll_expansion
avg_AIC_expansion = 2 * k - 2 * avg_ll_expansion
print "Average AIC:", avg_AIC_expansion
print datetime.now() - startTime
