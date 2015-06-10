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


# exponential growth model

growth_func = dadi.Demographics1D.growth
# params are nu: ratio of population size & T: time that change happened
growth_params = array([2,0.05])
growth_upper_bound = [100, 10]
growth_lower_bound = [1e-2, 0]
growth_func_ex = dadi.Numerics.make_extrap_log_func(growth_func)
growth_model = growth_func_ex(growth_params, ns, pts_l)
growth_ll = dadi.Inference.ll_multinom(growth_model, data)


growth_p0 = dadi.Misc.perturb_params(growth_params, fold=1,
                                        upper_bound = growth_upper_bound)

growth_popt = dadi.Inference.optimize_log(growth_p0, data, 
                                            growth_func_ex, pts_l,
                                            lower_bound = growth_lower_bound,
                                            upper_bound = growth_upper_bound,
                                            maxiter=3)
growth_model = growth_func_ex(growth_popt, ns, pts_l)
growth_ll_opt = dadi.Inference.ll_multinom(growth_model, data)

outfile.write("%f\t%f\t%f\n" % (growth_popt[0], growth_popt[1],
growth_ll_opt))
outfile.close()
