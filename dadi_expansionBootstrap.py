#!/usr/bin/env python

import sys
import numpy 
from numpy import array
import dadi


def usage():
    print "dadi_boostrapExpansion.py <input file> <output file>"

def likelihood_grid(function, data, ns, pts_l):
    outfile = open("likelihood_grid.txt", "w")
    outfile.write("nu\tT\tLL\n")
    for T in numpy.arange(0.001, 10, .1):
        for nu in numpy.arange(0.001, 10, .1):
            params = array([nu, T])
            model = function(params, ns, pts_l)
            ll = dadi.Inference.ll_multinom(model, data)
            outfile.write("%f\t%f\t%f\n" % (nu, T, ll))
    outfile.close()


if len(sys.argv) != 3:
    usage()
    sys.exit()

infile = sys.argv[1]
outfile = open(sys.argv[2], "w")

data = dadi.Spectrum.from_file(infile)
ns = data.sample_sizes

pts_l = [10,20,30] # grid point settings


# Instantaneous expansion model

expansion_func = dadi.Demographics1D.two_epoch
# params are nu: ratio of population size & T: time that change happened
expansion_params = array([2,0.05])
expansion_upper_bound = [100, 10]
expansion_lower_bound = [1e-2, 0]
expansion_func_ex = dadi.Numerics.make_extrap_log_func(expansion_func)
expansion_model = expansion_func_ex(expansion_params, ns, pts_l)
expansion_ll = dadi.Inference.ll_multinom(expansion_model, data)


expansion_p0 = dadi.Misc.perturb_params(expansion_params, fold=1,
                                        upper_bound = expansion_upper_bound)

expansion_popt = dadi.Inference.optimize_log(expansion_p0, data, 
                                            expansion_func_ex, pts_l,
                                            lower_bound = expansion_lower_bound,
                                            upper_bound = expansion_upper_bound,
                                            maxiter=3)
expansion_model = expansion_func_ex(expansion_popt, ns, pts_l)
expansion_ll_opt = dadi.Inference.ll_multinom(expansion_model, data)

outfile.write("%f\t%f\t%f\n" % (expansion_popt[0], expansion_popt[1],
expansion_ll_opt))
outfile.close()
