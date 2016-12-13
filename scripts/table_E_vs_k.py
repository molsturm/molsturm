#!/usr/bin/env python3

import common
import numpy as np
import sys
import os.path
import argparse

def add_args_to_parser(parser):
    """Add all required variables to the parser"""
    common.add_args_to_parser(parser)
    parser.add_argument("--kmin",metavar="kmin",type=float,default=None,help="Minimal value for k")
    parser.add_argument("--kmax",metavar="kmax", type=float,default=None,help="Maximal value for k")
    parser.add_argument("--count",metavar="n", type=int,default=100,help="Number of points to compute (default: 100)")
    parser.add_argument("--file",metavar="out",type=str,default=None,help="File to print computed data to (default: Z_nmax.dat)")

def parse_args(parser):
    """
    Parse all args from the parser object required for this module
    Return the args parsed, raise SystemExit on any error.
    """
    args = parser.parse_args()

    if args.kmin is None or args.kmax is None:
        raise SystemExit("Script needs --kmin, --kmax")

    if args.file is None:
        args.file=str(args.Z_charge) + "_" + str(args.n_max) + ".dat"

    if os.path.exists(args.file):
        raise SystemExit("File '"+args.file + "' exists already. Please remove or rename")

    return args

def energies(ks,params,print_progress=True):
    """ Compute energies for a set of k values
        and the common parameters params.

        if print_progress then a progress bar of the
        process is displayed.
    """
    if (print_progress):
        # setup a nice progress bar
        width=50
        sys.stdout.write("   ["+(" "*width)+"]")
        sys.stdout.flush()
        sys.stdout.write("\b" * (width+1))

        n_dash=0

    count = len(ks)

    Es=[ 0. for i in range(count)]
    for i in range(count):
        # Run molsturm with params:
        params["kexp"] = ks[i]
        Es[i] = common.molsturm(params)["hf_energy"]

        # Print progress:
        if print_progress and int(i*width/(count)) > n_dash:
            n_dash+=1
            sys.stdout.write("#")
            sys.stdout.flush()
    if print_progress:
        for i in range(n_dash,width):
            sys.stdout.write("#")
        sys.stdout.write("\n")
    return Es

#
# ---------------------------------------------------
#

if __name__ == "__main__":
    # Parse arguments:
    parser = argparse.ArgumentParser("Write a table of E versus k to a file")
    add_args_to_parser(parser)
    args = parse_args(parser)

    # Setup problem:
    params = common.get_params(args)
    ks = np.linspace(args.kmin,args.kmax,args.count+1)

    print("Running energy calculations:")
    Es = energies(ks,params,print_progress=True)

    print("Dumping data to", args.file)
    with open(args.file,"w") as f:
        f.write("# Z:       {0:6.2f}\n".format(args.Z_charge))
        f.write("# alpha:   {0:2d}\n".format(args.alpha))
        f.write("# beta:    {0:2d}\n".format(args.beta))
        f.write("#\n")
        f.write("# nmax:    {0:2d}\n".format(args.n_max))
        f.write("# lmax:    {0:2d}\n".format(args.l_max))
        f.write("#\n")
        f.write("# error:   {0:6.3E}\n".format(args.error))
        f.write("#\n")

        f.write("# {0:15s}    {1:15s}\n".format("k","E"))
        for i in range(len(ks)):
            f.write(" {0:15f}    {1:15f}\n".format(ks[i],Es[i]))

