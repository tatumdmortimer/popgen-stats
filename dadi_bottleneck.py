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

def rep_ll_bottleneck(data, ns, bottleneck_params, bottleneck_upper_bound,
                     bottleneck_lower_bound, bottleneck_func_ex, pts_l):
    bottleneck_ll_vals = []
    bottleneck_nuB_vals = []
    bottleneck_nuF_vals = []
    bottleneck_tauB_vals = []
    bottleneck_tauF_vals = []
    outfile = open("dadi_bottleneckModel.txt", "w")
    outfile.write("nuB\tnuF\ttauB\ttauF\tll\n")
    for i in range(100):
        bottleneck_p0 = dadi.Misc.perturb_params(bottleneck_params, fold=1,
                                        upper_bound = bottleneck_upper_bound)

        bottleneck_popt = dadi.Inference.optimize_log(bottleneck_p0, data,
                                            bottleneck_func_ex, pts_l,
                                            lower_bound = bottleneck_lower_bound,
                                            upper_bound = bottleneck_upper_bound,
                                            maxiter=100)
        bottleneck_model = bottleneck_func_ex(bottleneck_popt, ns, pts_l)
        bottleneck_ll_opt = dadi.Inference.ll_multinom(bottleneck_model, data)
        outfile.write(
            "{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                bottleneck_popt[0],
                bottleneck_popt[1],
                bottleneck_popt[2],
                bottleneck_popt[3],
                bottleneck_ll_opt))
        bottleneck_nuB_vals.append(bottleneck_popt[0])
        bottleneck_nuF_vals.append(bottleneck_popt[1])
        bottleneck_tauB_vals.append(bottleneck_popt[2])
        bottleneck_tauF_vals.append(bottleneck_popt[3])
        bottleneck_ll_vals.append(bottleneck_ll_opt)
    avg_ll_bottleneck = numpy.median(bottleneck_ll_vals)
    avg_nuB_bottleneck = numpy.median(bottleneck_nuB_vals)
    avg_nuF_bottleneck = numpy.median(bottleneck_nuF_vals)
    avg_tauB_bottleneck = numpy.median(bottleneck_tauB_vals)
    avg_tauF_bottleneck = numpy.median(bottleneck_tauF_vals)
    max_ll_index, max_ll_bottleneck = max(
        enumerate(bottleneck_ll_vals), key=lambda p: p[1])
    max_nuB_bottleneck = bottleneck_nuB_vals[max_ll_index]
    max_nuF_bottleneck = bottleneck_nuF_vals[max_ll_index]
    max_tauB_bottleneck = bottleneck_tauB_vals[max_ll_index]
    max_tauF_bottleneck = bottleneck_tauF_vals[max_ll_index]
    return (
        avg_ll_bottleneck,
        avg_nuB_bottleneck,
        avg_nuF_bottleneck,
        avg_tauB_bottleneck,
        avg_tauF_bottleneck,
        max_ll_bottleneck,
        max_nuB_bottleneck,
        max_nuF_bottleneck,
        max_tauB_bottleneck,
     max_tauF_bottleneck)
startTime = datetime.now()

pts_l = [110, 120, 130]
bottleneck_func = dadi.Demographics1D.three_epoch
# Params are nuB,nuF,TB,TF; nuB: Ratio of bottleneck population size to ancient pop size, nuF: Ratio of contemporary to ancient pop size,
# TB: Length of bottleneck and TF: Time since bottleneck recovery
bottleneck_params = array([2,2,0.05,0.05])
bottleneck_upper_bound = [100, 100, 10, 10]
bottleneck_lower_bound = [1e-2, 1e-2, 0, 0]
bottleneck_func_ex = dadi.Numerics.make_extrap_log_func(bottleneck_func)
k = len(bottleneck_params)

args = get_args()
data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()

avg_ll_bottleneck, avg_nuB_bottleneck, avg_nuF_bottleneck, avg_tauB_bottleneck, avg_tauF_bottleneck, max_ll_bottleneck, max_nuB_bottleneck, max_nuF_bottleneck, max_tauB_bottleneck, max_tauF_bottleneck = rep_ll_bottleneck(data, ns, bottleneck_params, bottleneck_upper_bound, bottleneck_lower_bound, bottleneck_func_ex, pts_l)

print "Average optimized log-likelihood:",avg_ll_bottleneck
avg_AIC_bottleneck = 2 * k - 2 * avg_ll_bottleneck
print "Average AIC:",avg_AIC_bottleneck

print "Maximum optimized log-likelihood:", max_ll_bottleneck
max_AIC_bottleneck = 2 * k - 2 * max_ll_bottleneck
print "Optimal AIC:", max_AIC_bottleneck

bottleneck_best_model = bottleneck_func_ex(
    [max_nuB_bottleneck, max_nuF_bottleneck, max_tauB_bottleneck, max_tauF_bottleneck], ns, pts_l)
bottleneck_sfs = dadi.Inference.optimally_scaled_sfs(
    bottleneck_best_model, data)

# Output SFS for bottlegrowth model
with open("bottleneckModelSFS.txt", "w") as sfsfile:
    for i in range(1, len(bottleneck_sfs)-1):
        sfsfile.write(str(bottleneck_sfs[i]) + "\n")
print datetime.now() - startTime
