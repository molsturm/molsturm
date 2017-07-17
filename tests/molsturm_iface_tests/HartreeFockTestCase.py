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

from NumCompTestCase import NumCompTestCase
import molsturm

class HartreeFockTestCase(NumCompTestCase):
  def compare_hf_0_energies(self, case, hfres):
    num_tol = case["testing"]["numeric_tolerance"]
    ref     = case["hf"]

    # Compute derived energies:
    hfene = molsturm.compute_derived_hartree_fock_energies(hfres)

    # Compare energies. Always compare energy_ground_state first
    for key in [ "energy_ground_state" ] + sorted(ref):
      if not key.startswith("energy_"):
        continue

      if key in hfres:
        computed = hfres[key]
      else:
        computed = hfene[key]

      self.assertAlmostEqual(computed, ref[key], tol=num_tol, prefix=key + ": ")


  def compare_hf_1_orben_orbcoeff(self, case, hfres):
    testing   = case["testing"]
    num_tol   = testing["numeric_tolerance"]
    orben_tol = testing.get("orben_tolerance", num_tol)
    ref       = case["hf"]

    for k in [ "orben_f", "orbcoeff_bf" ]:
      if k in ref:
        self.assertArrayAlmostEqual(hfres["orben_f"], ref["orben_f"], tol=orben_tol,
                                    prefix=k + ": ")
      else:
        print("Skipping HF test for " + k + " for " + case["testing"]["name"]
              + ", since no reference data available.")


  def compare_hf_2_convergence(self, case, hfres):
    testing          = case["testing"]
    params           = case["params"]
    num_tol          = testing["numeric_tolerance"]

    if "conv_tol" in params:
      self.assertLessEqual(hfres["final_error_norm"], params["conv_tol"])
    self.assertLessEqual(hfres["n_iter"], testing["max_n_iter"])


  def compare_hf_3_spin(self, case, hfres):
    testing          = case["testing"]
    num_tol          = testing["numeric_tolerance"]
    ref              = case["hf"]

    self.assertEqual(hfres["restricted"], ref["restricted"])
    spin_squared_tol = testing.get("spin_squared_tolerance", num_tol)
    self.assertAlmostEqual(hfres["spin_squared"], ref["spin_squared"], tol=spin_squared_tol,
                           prefix="spin squared: ")


  def compare_hf_results(self, case, hfres):
    """
    Call all compare functions after another """
    for fun in dir(HartreeFockTestCase):
      if fun == "compare_hf_results":
        continue # Do not call ourselves.

      if fun.startswith("compare_hf_"):
        getattr(self, fun) (case, hfres)

