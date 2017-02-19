#!/usr/bin/env python3

# TODO use chebychev polynomial interpolation with 4 or 6 points in a sensible vicinity
#      around a guess k and then optimise that polynomial instead.

import common
import math
from scipy.optimize import minimize
import argparse

def add_args_to_parser(parser):
    """Add all required variables to the parser"""
    common.add_args_to_parser(parser)
    parser.add_argument("--kerror",metavar="tol",type=float,default=1e-6,help="Maximal remaining error in k (Default: 1e-6)")
    parser.add_argument("--kguess",metavar="k", type=float,default=None,help="A guess for k (taken as sqrt(Z) if absent)")

def parse_args(parser):
    """
    Parse all args from the parser object required for this module
    Return the args parsed, raise SystemExit on any error.
    """
    args = parser.parse_args()

    if args.kguess is None:
        if args.Z_charge is None:
            args.kguess = 1.
        else:
            args.kguess = math.sqrt(args.Z_charge)

    return args

def find_minimal_k(params,kguess,kerror,disp=True):
    """
    Find the minimal k for the system and basis set
    described by param.
    Start from kguess and optimise until an error threshold
    below kerror is reached.

    disp:  Control whether to display convergence progress
    """
    # The energy functional to optimise
    def energy(args):
        params["kexp"] = args[0]
        return common.molsturm(params)["hf_energy"]

    if disp:
        # Setup printing of convergence:

        print("   iter    kexp")
        def priter(args):
            print("   {0:4d}    {1}".format(priter.count,args[0]))
            priter.count += 1
        priter.count=0
        priter([kguess])
    else:
        priter=None

    # Minimise the functional wrt. k
    res = minimize(energy,kguess,tol=args.kerror,options={"disp": disp},callback=priter)

    if not res.success:
        raise SystemExit("Error converging, got up to k="+str(res.x))

    return res.x[0]

#
# --------------------------------------------
#

if __name__ == "__main__":
    # Parse arguments:
    parser = argparse.ArgumentParser("Minimise the k exponent of a coulomb sturmian basis for a system.")
    add_args_to_parser(parser)
    args = parse_args(parser)
    params = common.get_params(args)

    kmin = find_minimal_k(params,args.kguess,args.kerror,disp=True)

    # Print results: n=Number of digits to print
    n = int(-math.log(args.kerror)/math.log(10)+0.5)
    print(("kopt = {0:."+str(n)+"f}  (tol={1})").format(kmin,args.kerror))
