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
import molsturm.posthf
import yaml
import numpy as np

def run_fci_for(atom, n_roots, experimental_ground_state, **kwargs):
  k_exp = np.sqrt(-experimental_ground_state*2)
  params = {
    "atoms":         atom,
    "basis_type":    "sturmian/atomic/cs_dummy",
    "n_max":         10,
    "l_max":         1,
    "k_exp":         k_exp,
    #
    "print_iterations": True,
    #
    "eigensolver":   "lapack",
    "guess_esolver": "lapack",
    "error":          1e-10,
    #
    "export_repulsion_integrals": True,
  }
  params.update(kwargs)

  #
  # Hartree-Fock
  #
  res = molsturm.hartree_fock(**params)
  molsturm.print_convergence_summary(res)
  molsturm.print_energies(res)
  molsturm.print_mo_occupation(res)

  #
  # Full-CI
  #
  print("\nRunning FCI ... please wait\n")
  res_fci = molsturm.posthf.fci(res, verbosity=5, n_roots=n_roots)

  #
  # Print results
  #
  def spin_string(multiplity):
    m = { 1: "S", 2: "D", 3: "T" }
    for k in m:
      if abs(float(multiplity) - float(k)) < 1e-6:
        return m[k]
    else:
      raise ValueError("Unknown multiplity: " + str(multiplity))

  print("\nFCI energies:")
  E_GS = res_fci["states"][0]["energy"]
  eV =  27.21138602  # au zu eV
  count_system=dict()
  for i in range(n_roots):
    spin = spin_string(res_fci["states"][i]["multiplicity"])
    cnt = count_system.get(spin,0)
    count_system[spin] = cnt+1
    en = res_fci["states"][i]["energy"]
    print( "{0:2d} {1:3s}  {2:7.5g} a.u.   Δ = {3:7.5g} a.u. = {4:7.5g} eV".format(
      i, spin + str(cnt), en, en - E_GS, eV*(en - E_GS)))

  with open(atom + "_cs_fci.yaml", "w") as f:
    yaml.safe_dump([{
        "energy":        float(res_fci["states"][i]["energy"]),
        "multiplicity":  float(res_fci["states"][i]["multiplicity"]),
      } for i in range(n_roots) ], f)

  print()

  print("Experiment GS:", experimental_ground_state)
  print("          ΔGS:", E_GS - experimental_ground_state)

  molsturm.print_quote(res)

