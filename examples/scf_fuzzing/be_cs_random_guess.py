#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the molsturm authors
##
## This file is part of molsturm.
##
## molsturm is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## molsturm is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with molsturm. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

import gint
import molsturm
import multiprocessing
import numpy as np

# This scripts runs `n_points` molsturm calculations starting from
# a random guess and checks the SCF energies we get if we converge
# up to an error of error
n_points = 100
error = 1e-7

system = molsturm.MolecularSystem("Be")
basis = gint.sturmian.atomic.Basis(system, k_exp=1.3, n_max=3, l_max=2,
                                   backend="cs_static14")

# We want the data gathering to be parallel ...
try:
    nproc = multiprocessing.cpu_count()
except NotImplementedError:
    nproc = 2  # arbitrary default


# Function which runs molsturm and sends the total energy
# and the number of iterations to the main process
def run_molsturm(x):
    res = molsturm.hartree_fock(system, basis, eigensolver="lapack", max_iter=200,
                                conv_tol=error, guess="random")
    return (res["n_iter"], res["energy_ground_state"])


# Use tqdm to keep track of the progress if it is installed.
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it, **kwargs):
        return it

# Run SCF in parallel and extract the values we care about:
print("Running SCF with", n_points, "random guesses on ", nproc,
      " processes.", flush=True)
with multiprocessing.Pool(processes=nproc) as pool:
    results = [
        (res[0], res[1], np.round(res[1], decimals=-int(np.log10(error))))
        for res in tqdm(pool.imap_unordered(run_molsturm, range(n_points)),
                        total=n_points)
    ]

# Extract raw energy values and bin the results we obtained.
energies = np.array([res[1] for res in results])
binning = dict()
for r in results:
    element = (r[0], r[2])
    try:
        binning[element] += 1
    except KeyError:
        binning[element] = 1


fmt = "    {0:5}  {1:5}  {2}"
print(fmt.format("count", "niter", "energy"))
for b in sorted(binning, key=lambda x: x[1]):
    print(fmt.format(binning[b], b[0], b[1]))
print()
print("Mean energy:      ", np.average(energies), "  stddev ",
      np.sqrt(np.var(energies)))
