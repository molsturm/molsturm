import subprocess
import re

binary="../build/examples/hf/hf"

def getE(params):
    res = subprocess.check_output(
        [binary,"--Z", str(params["Z"]), 
         "--alpha", str(params["alpha"]), 
         "--beta", str(params["beta"]), 
         "--nmax", str(params["nmax"]), 
         "--kexp", str(params["k"])]
    ).splitlines()

    expr = re.compile("E_total *= *(-?[0-9.eE]+)")
    for line in res:
        foundit = expr.search(str(line))
        if foundit:
            return float(foundit.group(1))
    raise SystemExit("something wrong here")

