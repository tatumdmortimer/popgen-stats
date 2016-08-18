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

startTime = datetime.now()

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Run dadi analysis')
    parser.add_argument("sfs", help="SFS in dadi format",
                    type=argparse.FileType('r'))
    parser.add_argument("-f", "--fold", help="Fold the SFS before analysis",
                        action="store_true")
    return parser.parse_args()

def rep_ll_bottlegrowth(data, ns, bottlegrowth_params, bottlegrowth_upper_bound,
                     bottlegrowth_lower_bound, bottlegrowth_func_ex, pts_l):
    bottlegrowth_ll_vals = []
    bottlegrowth_nuB_vals = []
    bottlegrowth_nuF_vals = []
    bottlegrowth_tau_vals = []
    outfile = open("dadi_bottlegrowthModel.txt", "w")
    outfile.write("nuB\tnuF\ttau\tll\n")
    for i in range(100):
        bottlegrowth_p0 = dadi.Misc.perturb_params(bottlegrowth_params, fold=1,
                                        upper_bound = bottlegrowth_upper_bound)

        bottlegrowth_popt = dadi.Inference.optimize_log(bottlegrowth_p0, data,
                                            bottlegrowth_func_ex, pts_l,
                                            lower_bound = bottlegrowth_lower_bound,
                                            upper_bound = bottlegrowth_upper_bound,
                                            maxiter=100)
        bottlegrowth_model = bottlegrowth_func_ex(bottlegrowth_popt, ns, pts_l)
        bottlegrowth_ll_opt = dadi.Inference.ll_multinom(bottlegrowth_model, data)
        outfile.write(
            "{0}\t{1}\t{2}\t{3}\n".format(
                bottlegrowth_popt[0],
                bottlegrowth_popt[1],
                bottlegrowth_popt[2],
                bottlegrowth_ll_opt))
        bottlegrowth_nuB_vals.append(bottlegrowth_popt[0])
        bottlegrowth_nuF_vals.append(bottlegrowth_popt[1])
        bottlegrowth_tau_vals.append(bottlegrowth_popt[2])
        bottlegrowth_ll_vals.append(bottlegrowth_ll_opt)
    avg_ll_bottlegrowth = numpy.median(bottlegrowth_ll_vals)
    avg_nuB_bottlegrowth = numpy.median(bottlegrowh_nuB_vals)
    avg_nuF_bottlegrowth = numpy.median(bottlegrowh_nuF_vals)
    avg_tau_bottlegrowth = numpy.median(expansion_tau_vals)
    max_ll_index, max_ll_bottlegrowth = max(
        enumerate(bottlegrowth_ll_vals), key=lambda p: p[1])
    max_nuB_bottlegrowth = bottlegrowth_nuB_vals[max_ll_index]
    max_nuF_bottlegrowth = bottlegrowth_nuB_vals[max_ll_index]
    max_tau_bottlegrowth = bottlegrowth_tau_vals[max_ll_index]
    return (
        avg_ll_bottlegrowth,
        avg_nuB_bottlegrowth,
        avg_nuF_bottlegrowth,
        avg_tau_bottlegrowth,
        max_ll_bottlegrowth,
        max_nuB_bottlegrowth,
        max_nuF_bottlegrowth,
     max_tau_bottlegrowth)

startTime = datetime.now()

pts_l = [110, 120, 130]
bottlegrowth_func = dadi.Demographics1D.bottlegrowth
# Params are nuB,nuF,T
bottlegrowth_params = array([2,2,0.05])
bottlegrowth_upper_bound = [100, 100, 10]
bottlegrowth_lower_bound = [1e-2, 1e-2, 0]
bottlegrowth_func_ex = dadi.Numerics.make_extrap_log_func(bottlegrowth_func)
bottlegrowth_model = bottlegrowth_func_ex(bottlegrowth_params, ns, pts_l)
k = len(bottlegrowth_params)

args = get_args()
data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()

avg_ll_bottlegrowth, avg_nuB_bottlegrowth, avg_nuF_bottlegrowth, avg_tau_bottlegrowth, max_ll_bottlegrowth, max_nuB_bottlegrowth, max_nuF_bottlegrowth, max_tau_bottlegrowth =
    rep_ll_bottlegrowth(data, ns, bottlegrowth_params, bottlegrowth_upper_bound, bottlegrowth_lower_bound, bottlegrowth_func_ex, pts_l)

print "Average optimized log-likelihood:",avg_ll_bottlegrowth
avg_AIC_bottlegrowth = 2 * k - 2 * avg_ll_bottlegrowth
print "Average AIC:",avg_AIC_bottlegrowth

print "Maximum optimized log-likelihood:", max_ll_bottlegrowth
max_AIC_bottlegrowth = 2 * k - 2 * max_ll_bottlegrowth
print "Optimal AIC:", max_AIC_bottlegrowth

bottlegrowth_best_model = bottlegrowth_func_ex(
    [max_nuB_bottlegrowth, max_nuF_bottlegrowth, max_tau_bottlegrowth], ns, pts_l)
bottlegrowth_sfs = dadi.Inference.optimally_scaled_sfs(
    bottlegrowth_best_model, data)

# Output SFS for bottlegrowth model
with open("bottlegrowthModelSFS.txt", "w") as sfsfile:
    for i in range(1, len(bottlegrowth_sfs)-1):
        sfsfile.write(str(bottlegrowth_sfs[i]) + "\n")
print datetime.now() - startTime
