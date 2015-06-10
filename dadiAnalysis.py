#!/usr/bin/env python

import sys
import argparse
from Bio.Phylo.PAML.chi2 import cdf_chi2
import numpy 
from numpy import array
import dadi

def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Run dadi analysis')
    parser.add_argument("sfs", help="SFS in dadi format",
                    type=argparse.FileType('r'))
    parser.add_argument("-f", "--fold", help="Fold the SFS before analysis",
                        action="store_true")
    return parser.parse_args()

def likelihood_grid(function, data, ns, pts_l, func_name):
    outfile = open("likelihood_grid_{0}.txt".format(func_name), "w")
    outfile.write("nu\tT\tLL\n")
    for T in numpy.arange(0.001, 1, 0.01):
        for nu in numpy.arange(0.001, 50, 1):
            params = array([nu, T])
            model = function(params, ns, pts_l)
            ll = dadi.Inference.ll_multinom(model, data)
            outfile.write("%f\t%f\t%f\n" % (nu, T, ll))
    outfile.close()

args = get_args()

data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()
print "Number of samples: %s" % ns

thetaW = data.Watterson_theta()
print "Watterson's theta: %f" % thetaW

pi = data.pi()
print "Pi: %f" % pi

D = data.Tajima_D()
print "Tajima's D: %f" % D

pts_l = [10,20,30] # grid point settings


# Neutral model

neutral_func = dadi.Demographics1D.snm
neutral_params = array([])
neutral_upper_bound = []
neutral_func_ex = dadi.Numerics.make_extrap_log_func(neutral_func)
neutral_model = neutral_func_ex(neutral_params, ns, pts_l)
neutral_ll = dadi.Inference.ll_multinom(neutral_model, data)

print "Neutral model log-likelihood: %f" % neutral_ll

# Instantaneous expansion model

expansion_func = dadi.Demographics1D.two_epoch
# params are nu: ratio of population size & T: time that change happened
expansion_params = array([2,0.05])
expansion_upper_bound = [100, 10]
expansion_lower_bound = [1e-2, 0]
expansion_func_ex = dadi.Numerics.make_extrap_log_func(expansion_func)
expansion_model = expansion_func_ex(expansion_params, ns, pts_l)
expansion_ll = dadi.Inference.ll_multinom(expansion_model, data)

print "Expansion model log-likelihood: %f" % expansion_ll


expansion_p0 = dadi.Misc.perturb_params(expansion_params, fold=1,
                                        upper_bound = expansion_upper_bound)

expansion_popt = dadi.Inference.optimize_log(expansion_p0, data, 
                                            expansion_func_ex, pts_l,
                                            lower_bound = expansion_lower_bound,
                                            upper_bound = expansion_upper_bound,
                                            maxiter=3)
print "Optimized parameters", repr(expansion_popt)
expansion_model = expansion_func_ex(expansion_popt, ns, pts_l)
expansion_ll_opt = dadi.Inference.ll_multinom(expansion_model, data)
print "Optimized log-likelihood:", expansion_ll_opt

# Exponential growth model

growth_func = dadi.Demographics1D.growth
# params are nu: ratio of population size & T: time that change happened
growth_params = array([2,0.05])
growth_upper_bound = [100, 10]
growth_lower_bound = [1e-2, 0]
growth_func_ex = dadi.Numerics.make_extrap_log_func(growth_func)
growth_model = growth_func_ex(growth_params, ns, pts_l)
growth_ll = dadi.Inference.ll_multinom(growth_model, data)

print "Exponential growth model log-likelihood: %f" % growth_ll

growth_p0 = dadi.Misc.perturb_params(growth_params, fold=1,
                                        upper_bound = growth_upper_bound)

growth_popt = dadi.Inference.optimize_log(growth_p0, data, 
                                            growth_func_ex, pts_l,
                                            lower_bound = growth_lower_bound,
                                            upper_bound = growth_upper_bound,
                                            maxiter=3)
print "Optimized parameters", repr(growth_popt)
growth_model = growth_func_ex(growth_popt, ns, pts_l)
growth_ll_opt = dadi.Inference.ll_multinom(growth_model, data)
print "Optimized log-likelihood:", growth_ll_opt

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

# Output SFS for expansion model
expansion_sfs = dadi.Inference.optimally_scaled_sfs(expansion_model, data)
expansion_sfs_file = open("expansionModelSFS.txt", 'w')
for i in range(1,len(expansion_sfs)-1):
    expansion_sfs_file.write(str(expansion_sfs[i]) + '\n')
expansion_sfs_file.close()

# Output SFS for growth model
growth_sfs = dadi.Inference.optimally_scaled_sfs(growth_model, data)
growth_sfs_file = open("growthModelSFS.txt", 'w')
for i in range(1,len(growth_sfs)-1):
    growth_sfs_file.write(str(growth_sfs[i]) + '\n')
growth_sfs_file.close()

if expansion_ll_opt > growth_ll_opt:
    print "Testing significance of expansion..."
    LRTstat_neutral = 2*(expansion_ll_opt - neutral_ll)
    degrees = len(expansion_params)
    print "LRT Statistic- Neutral:", LRTstat_neutral
    p_neutral = cdf_chi2(degrees, LRTstat_neutral)
    print "p-value neutral =",p_neutral
    LRTstat_growth = 2*(expansion_ll_opt - growth_ll_opt)
    print "LRT Statistic- Exponential Growth:", LRTstat_growth
    p_growth = cdf_chi2(degrees, LRTstat_growth)
    print "p-value growth =",p_growth
    if p_neutral < 0.05:
        print "Working on likelihood surface..."
        likelihood_grid(expansion_func_ex, data, ns, pts_l, "expansion")
        if p_growth > 0.05:
            likelihood_grid(growth_func_ex, data, ns, pts_l, "growth")

else:
    print "Testing significance of exponential growth..."
    LRTstat_neutral = 2*(growth_ll_opt - neutral_ll)
    degrees = len(growth_params)
    print "LRT Statistic - Neutral:", LRTstat_neutral
    p_neutral = cdf_chi2(degrees, LRTstat_neutral)
    print "p-value neutral=",p_neutral
    LRTstat_expansion = 2*(growth_ll_opt - expansion_ll_opt)
    print "LRT Statistic - Expansion:", LRTstat_expansion
    p_expansion = cdf_chi2(degrees, LRTstat_expansion)
    print "p-value expansion =", p_expansion
    if p_neutral < 0.05:
        print "Working on likelihood surface..."
        likelihood_grid(growth_func_ex, data, ns, pts_l, "growth")
        if p_expansion > 0.05:
            likelihood_grid(expansion_func_ex, data, ns, pts_l, "expansion")

