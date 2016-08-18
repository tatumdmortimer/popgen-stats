#!/usr/bin/env python

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


def rep_ll_expansion(data, ns, expansion_params, expansion_upper_bound,
                     expansion_lower_bound, expansion_func_ex, pts_l):
    expansion_ll_vals = []
    expansion_nu_vals = []
    expansion_tau_vals = []
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
        expansion_ll_vals.append(expansion_ll_opt)
        expansion_nu_vals.append(expansion_popt[0])
        expansion_tau_vals.append(expansion_popt[1])
    outfile.close()
    avg_ll_expansion = numpy.median(expansion_ll_vals)
    avg_nu_expansion = numpy.median(expansion_nu_vals)
    avg_tau_expansion = numpy.median(expansion_tau_vals)
    max_ll_index, max_ll_expansion = max(
        enumerate(expansion_ll_vals), key=lambda p: p[1])
    max_nu_expansion = expansion_nu_vals[max_ll_index]
    max_tau_expansion = expansion_tau_vals[max_ll_index]
    return (
        avg_ll_expansion,
        avg_nu_expansion,
        avg_tau_expansion,
        max_ll_expansion,
        max_nu_expansion,
     max_tau_expansion)


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

avg_ll_expansion, avg_nu_expansion, avg_tau_expansion, max_ll_expansion, max_nu_expansion, max_tau_expansion = rep_ll_expansion(
     data, ns, expansion_params, expansion_upper_bound, expansion_lower_bound, expansion_func_ex, pts_l)
print "Average optimized log-likelihood:", avg_ll_expansion
avg_AIC_expansion = 2 * k - 2 * avg_ll_expansion
print "Average AIC:", avg_AIC_expansion

print "Maximum optimized log-likelihood:", max_ll_expansion
max_AIC_expansion = 2 * k - 2 * max_ll_expansion
print "Optimal AIC:", max_AIC_expansion

expansion_best_model = expansion_func_ex(
    [max_nu_expansion, max_tau_expansion], ns, pts_l)
expansion_sfs = dadi.Inference.optimally_scaled_sfs(
    expansion_best_model, data)
with open("expansionModelSFS.txt", "w") as sfsfile:
    for i in range(1, len(expansion_sfs)-1):
        sfsfile.write(str(expansion_sfs[i]) + "\n")
print datetime.now() - startTime
