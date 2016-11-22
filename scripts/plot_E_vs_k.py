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
deltak=0.01
count=200

print("#  k       E")
for i in range(count):
    k = kmin+deltak*i
    params["k"] = k
    E = common.getE(params)
    print(k, "    ", E)

