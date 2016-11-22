#!/usr/bin/env python3

import common

params = {
    "Z": 4, 
    "alpha": 2,
    "beta": 2,
    "nmax": 3,
}

# search space
kmin=0.5
kmax=2.4

tolk = 1e-8

print("This is not the best algorithm and does not always converge to the actual minimum")

#
# --------------------------------------------
#

def getE(k):
    print("Running with k =",k)
    params["k"] = k
    E = common.getE(params)
    print("    E =", E)

Ekmin = getE(kmin)
Ekmax = getE(kmax)

while (abs(kmin-kmax) > tolk):
    kcurr= (kmin+kmax)/2
    Ecurr = getE(kcurr)
    if Ecurr < Ekmin:
        Ekmin = Ecurr
        kmin = kcurr
    elif Ecurr < Ekmax:
        Ekmax = Ecurr
        kmax = kcurr
    else:
        raise SystemExit("something wrong")

kcurr= (kmin+kmax)/2
Ecurr = getE(kcurr)
print("k: ", kcurr,"    E:", Ecurr)
