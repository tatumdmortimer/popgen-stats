#!/usr/bin/env python

import argparse
import numpy
from numpy import array
import dadi
import sys


def get_args():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Run dadi analysis')
    parser.add_argument("sfs", help="SFS in dadi format",
                        type=argparse.FileType('r'))
    parser.add_argument(
        "model",
        help="Demographic model to calculate likelihood for",
        choices=[
            'expansion',
            'growth',
            'bottleneck',
         'bottlegrowth'])
    parser.add_argument("-f", "--fold", help="Fold the SFS before analysis",
                        action="store_true")
    parser.add_argument(
        "-n",
        "--min",
        help="Minimum values for each parameter in surface",
     nargs='+', type=float)
    parser.add_argument(
        "-x",
        "--max",
        help="Maximum values for each parameter in surface",
     nargs='+', type=float)
    return parser.parse_args()


def likelihood_grid(function, data, ns, pts_l, func_name, n, x, s):
    outfile = open("likelihood_grid_{0}.txt".format(func_name), "w")
    outfile.write("nu\tT\tLL\n")
    for T in numpy.arange(n[0], x[0], s[0]):
        for nu in numpy.arange(n[1], x[1], s[1]):
            params = array([nu, T])
            model = function(params, ns, pts_l)
            ll = dadi.Inference.ll_multinom(model, data)
            outfile.write("%f\t%f\t%f\n" % (nu, T, ll))
    outfile.close()


def likelihood_grid_bottleneck(function, data, ns, pts_l, func_name, n, x, s):
    outfile = open("likelihood_grid_{0}.txt".format(func_name), "w")
    outfile.write("nuB\tnuF\tTB\tTF\tLL\n")
    for TF in numpy.arange(n[0], x[0], s[0]):
        for TB in numpy.arange(n[1], x[1], s[1]):
            for nuF in numpy.arange(n[2], x[2], s[2]):
                for nuB in numpy.arange(n[3], x[3], s[3]):
                    params = array([nuB, nuF, TB, TF])
                    model = function(params, ns, pts_l)
                    ll = dadi.Inference.ll_multinom(model, data)
                    outfile.write(
                        "%f\t%f\t%f\t%f\t%f\n" %
                        (nuB, nuF, TB, TF, ll))
    outfile.close()


def likelihood_grid_bottlegrowth(function, data, ns, pts_l, func_name):
    outfile = open("likelihood_grid_{0}.txt".format(func_name), "w")
    outfile.write("nuB\tnuF\tTF\tLL\n")
    for TF in numpy.arange(n[0], x[0], s[0]):
        for nuF in numpy.arange(n[1], x[1], s[1]):
            for nuB in numpy.arange(n[2], x[2], s[2]):
                params = array([nuB, nuF, TF])
                model = function(params, ns, pts_l)
                ll = dadi.Inference.ll_multinom(model, data)
                outfile.write("%f\t%f\t%f\t%f\n" % (nuB, nuF, TF, ll))
    outfile.close()


args = get_args()

data = dadi.Spectrum.from_file(args.sfs)
ns = data.sample_sizes
if args.fold:
    data = data.fold()
pts_l = [110, 120, 130]  # grid point settings

# check that bounds for likelihood surface make sense
if len(args.min) != len(args.max):
    print("Your minimum bounds and maximum bounds arguments are not the same length")
    sys.exit(1)
for i, n in enumerate(args.min):
    if n > args.max[i]:
        print(
            "The minimum bound for the likelihood surface is greater than the maximum bound for argument {0}".format(
                i +
                1))
        sys.exit(1)

# create step size for each parameter
step = []
for j, x in enumerate(args.max):
    r = x - args.min[j]
    step.append(r/50)
if args.model == "expansion":
    if len(args.min) != 2 or len(args.max) != 2:
        print("You did not give the correct number of parameters for the minimum or maximum argument for this model (2)")
        sys.exit(1)
    expansion_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)
    likelihood_grid(
        expansion_func_ex,
        data,
        ns,
        pts_l,
        "expansion",
        args.min,
        args.max,
     step)
elif args.model == "growth":
    if len(args.min) != 2 or len(args.max) != 2:
        print("You did not give the correct number of parameters for the minimum or maximum argument for this model (2)")
        sys.exit(1)
    growth_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.growth)
    likelihood_grid(
        growth_func_ex,
        data,
        ns,
        pts_l,
        "growth",
        args.min,
        args.max,
     step)
elif args.model == "bottleneck":
    if len(args.min) != 4 or len(args.max) != 4:
        print("You did not give the correct number of parameters for the minimum or maximum argument for this model (4)")
        sys.exit(1)
    bottleneck_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.three_epoch)
    likelihood_grid_bottleneck(
        bottleneck_func_ex,
        data,
        ns,
        pts_l,
     "bottleneck", args.min, args.max, step)
elif args.model == "bottlegrowth":
    if len(args.min) != 3 or len(args.max) != 3:
        print("You did not give the correct number of parameters for the minimum or maximum argument for this model (3)")
        sys.exit(1)
    bottlegrowth_func_ex = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.bottelgrowth)
    likelihood_grid_bottlegrowth(
        bottlegrowth_func_ex,
        data,
        ns,
        pts_l,
     "bottlegrowth", args.min, args.max, step)
