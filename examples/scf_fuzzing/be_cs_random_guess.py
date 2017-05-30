#!/usr/bin/env python3
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
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab

import molsturm
import multiprocessing
import numpy as np

# This scripts runs `n_points` molsturm calculations starting from
# a random guess and checks the SCF energies we get if we converge
# up to an error of error
n_points = 100
error = 1e-7

# Parameters for the molsturm caluclation:
params = {
  "atom_numbers": [4],
  "coords":       [[0,0,0]],
  #
  "basis_type":   "atomic/cs_dummy",
  "n_max":        3,
  "l_max":        2,
  "k_exp":        1.3,
  #
  "eigensolver":   "lapack",
  "guess_esolver": "lapack",
  #
  "max_iter":     100,
  "error":        error,
}
# ----------------------------------------------------

# We want the data gathering to be parallel ...
try:
  nproc = multiprocessing.cpu_count()
except NotImplementedError:
  nproc = 2 # arbitrary default

# Function which runs molsturm and sends the total energy
# and the number of iterations to the main process
def run_molsturm(x):
  amended = dict(params)
  amended["guess_method"] = "random"
  res = molsturm.hartree_fock(**amended)
  return (res["energy_total"], res["n_iter"])

# Use tqdm to keep track of the progress if it is installed.
try:
  from tqdm  import tqdm
except ImportError:
  def tqdm(it):
    return it

# Run SCF in parallel and extract the values we care about:
print("Running SCF with",n_points,"random guesses on ",nproc," processes.",flush=True)
with multiprocessing.Pool(processes=nproc) as pool:
  results = tqdm(pool.imap_unordered(run_molsturm,range(n_points)),total=n_points)

  # Separate the data by kind:
  results = [ i for i in zip(*results) ]
  energy = np.array(results[0])
  n_iter = np.array(results[1])

# Bin the energy into unique values under the error we used above
unique_values = np.unique(energy.round(decimals=-int(np.log10(error))))
binning = [ (val, np.sum( abs(energy - val) < error ))
            for val in unique_values ]

print()
print("Mean energy:      ",np.average(energy), "  stddev ",
      np.sqrt(np.var(energy)))
print("Mean number of iterations:   ", np.average(n_iter), " stddev ",
      np.sqrt(np.var(n_iter)))
print()
print("Unique energy values:")
for b in binning:
  print("    ",b)
