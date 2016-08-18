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


def rep_ll_growth(data, ns, growth_params, growth_upper_bound,
                     growth_lower_bound, growth_func_ex, pts_l):
    growth_ll_vals = []
    growth_nu_vals = []
    growth_tau_vals = []
    outfile = open("dadi_growthModel.txt", "w")
    outfile.write("nu\ttau\tll\n")
    for i in range(100):
        growth_p0 = dadi.Misc.perturb_params(
            growth_params, fold=1, upper_bound=growth_upper_bound)

        growth_popt = dadi.Inference.optimize_log(
            growth_p0, data, growth_func_ex, pts_l,
            lower_bound=growth_lower_bound,
            upper_bound=growth_upper_bound, maxiter=100)
        growth_model = growth_func_ex(growth_popt, ns, pts_l)
        growth_ll_opt = dadi.Inference.ll_multinom(growth_model, data)
        outfile.write(
            "{0}\t{1}\t{2}\n".format(
                growth_popt[0],
                growth_popt[1],
                growth_ll_opt))
        growth_ll_vals.append(growth_ll_opt)
        growth_nu_vals.append(growth_popt[0])
        growth_tau_vals.append(growth_popt[1])
    outfile.close()
    avg_ll_growth = numpy.median(growth_ll_vals)
    avg_nu_growth = numpy.median(growth_nu_vals)
    avg_tau_growth = numpy.median(growth_tau_vals)
    max_ll_index, max_ll_growth = max(
        enumerate(growth_ll_vals), key=lambda p: p[1])
    max_nu_growth = growth_nu_vals[max_ll_index]
    max_tau_growth = growth_tau_vals[max_ll_index]
    return (
        avg_ll_growth,
        avg_nu_growth,
        avg_tau_growth,
        max_ll_growth,
        max_nu_growth,
     max_tau_growth)


startTime = datetime.now()

pts_l = [110, 120, 130]
growth_func = dadi.Demographics1D.growth
growth_params = array([2, 0.05])
growth_upper_bound = [100, 10]
growth_lower_bound = [1e-2, 0]
growth_func_ex = dadi.Numerics.make_extrap_log_func(growth_func)
k = len(growth_params)

args = get_args()
data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()

avg_ll_growth, avg_nu_growth, avg_tau_growth, max_ll_growth, max_nu_growth, max_tau_growth = rep_ll_growth(
     data, ns, growth_params, growth_upper_bound, growth_lower_bound, growth_func_ex, pts_l)
print "Average optimized log-likelihood:", avg_ll_growth
avg_AIC_growth = 2 * k - 2 * avg_ll_growth
print "Average AIC:", avg_AIC_growth

print "Maximum optimized log-likelihood:", max_ll_growth
max_AIC_growth = 2 * k - 2 * max_ll_growth
print "Optimal AIC:", max_AIC_growth

growth_best_model = growth_func_ex(
    [max_nu_growth, max_tau_growth], ns, pts_l)
growth_sfs = dadi.Inference.optimally_scaled_sfs(
    growth_best_model, data)
with open("growthModelSFS.txt", "w") as sfsfile:
    for i in range(1, len(growth_sfs)-1):
        sfsfile.write(str(growth_sfs[i]) + "\n")
print datetime.now() - startTime
