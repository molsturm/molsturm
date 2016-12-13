import subprocess
import re

binary="../build/examples/hf/hf"

def molsturm(params):
    plist = [binary]
    for (k,v) in params.items():
        plist.append("--" + k)
        plist.append(str(v))

    res = subprocess.check_output(plist).splitlines()

    # Result dictionary
    out = {}

    expr = re.compile("E_total *= *(-?[0-9.eE]+)")
    for line in res:
        foundit = expr.search(str(line))
        if foundit:
            out["hf_energy"] = float(foundit.group(1))
            return out
    raise SystemExit("something wrong here")

def add_args_to_parser(parser):
    """
    Add the arguments understood by molsturm to the parser
    """
    parser.add_argument("--Z_charge",metavar="Z",type=float,default=None,help="Nuclear charge")
    parser.add_argument("--alpha",metavar="na",type=int,default=None,help="Number of alpha electrons")
    parser.add_argument("--beta",metavar="nb",type=int,default=None,help="Number of beta electrons")

    parser.add_argument("--error",metavar="tol", type=float,default=1e-4,help="Error for the computation (Default: 1e-4)")
    parser.add_argument("--max_iter",metavar="iter", type=float,default=100,help="Maximal number of iterations (Default: 100)")

    parser.add_argument("--n_max",metavar="n", type=int,default=3,help="Maximal number for n")
    parser.add_argument("--l_max",metavar="l", type=int,default=None,help="Maximal number for l (Default: n-1)")

def get_params(args):
    """
    Get the params for system, basis and convergence
    """
    params = {}

    # Basis Set params
    params["basis_type"] = "cs_naive"
    params["n_max"] = args.n_max

    if args.l_max is None:
        params["l_max"] = args.n_max -1
        args.l_max = args.n_max -1
    else:
        params["l_max"] = args.l_max

    if args.Z_charge is None or args.alpha is None or args.beta is None:
        raise SystemExit("Script needs --Z_charge, --alpha, --beta")

    # System params
    params["Z_charge"] = args.Z_charge
    params["alpha"] = args.alpha
    params["beta"] = args.beta

    # Convergence params
    params["n_eigenpairs"] = 100  #TODO To make molsturm converge
    params["error"] = args.error
    params["max_iter"] = args.max_iter

    return params

