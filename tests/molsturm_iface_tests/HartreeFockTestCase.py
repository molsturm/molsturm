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
from molsturm._basis import has_real_harmonics
from molsturm import INPUT_PARAMETER_KEY
import molsturm
import numpy as np
import time
import test_constants

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
    self.assertAlmostEqual(hfres["spin_squared"], ref["spin_squared"],
                           tol=spin_squared_tol,
                           prefix="spin squared: ")


  def compare_hf_4_hcore_symmetry(self, case, hfres):
    num_tol  = case["testing"]["numeric_tolerance"]
    num_zero = max(1e-14, num_tol/100)

    if not "hcore_ff" in hfres:
      print("Skipping hcore symmetry test for " + case["testing"]["name"]
            + ", since hcore_ff not available")
      return

    hcore = hfres["hcore_ff"]

    # Permutational symmetry
    diff = np.transpose(hcore) - hcore
    max_asym = np.max(np.abs(diff))
    self.assertLessEqual(max_asym, num_zero)

    # Spin symmetry
    noa = hfres["n_orbs_alpha"]
    max_ab = np.max(np.abs(hcore[:noa,noa:]))
    max_ba = np.max(np.abs(hcore[noa:,:noa]))
    self.assertLessEqual(max_ab, num_zero)
    self.assertLessEqual(max_ba, num_zero)


  def compare_hf_5_hcore_e1e(self, case, hfres):
    num_tol  = case["testing"]["numeric_tolerance"]

    if not "hcore_ff" in hfres:
      print("Skipping hcore e1e test for " + case["testing"]["name"]
            + ", since hcore_ff not available")
      return

    # Sum diagonal elements of hcore which make up the one-electron energy
    hcore = hfres["hcore_ff"]
    noa   = hfres["n_orbs_alpha"]
    na    = hfres["n_alpha"]
    nb    = hfres["n_beta"]
    comp_e1e = np.sum(hcore.diagonal()[:na]) + np.sum(hcore.diagonal()[noa:noa + nb])

    # Check against the computed one electron energy
    hfene = molsturm.compute_derived_hartree_fock_energies(hfres)
    self.assertAlmostEqual(hfene["energy_1e"], comp_e1e, tol=num_tol,
                           prefix="energy_1e from hcore: ")


  def compare_hf_6_eri_symmetry(self, case, hfres):
    num_tol  = case["testing"]["numeric_tolerance"]
    num_zero = max(1e-14, num_tol/100)

    if not "eri_ffff" in hfres:
      print("Skipping eri symmetry test for " + case["testing"]["name"]
            + ", since eri_ffff not available")
      return

    #
    # Test permutational symmetry
    #
    jirep = hfres["eri_ffff"] # in Mulliken notation
    allowed_permutations=[
      (1,0,3,2),  # (1) Swap indices in both shell pairs (TODO only for real)
      (2,3,0,1),  # (2) Swap shell pairs
      (3,2,1,0),  # Both (1) and (2)
    ]

    if has_real_harmonics(**hfres[INPUT_PARAMETER_KEY]):
      allowed_permutations.extend([
        (1,0,2,3),  # (3) Permutation inside 1st shell pair
        (0,1,3,2),  # (4) Permutation inside 2nd shell pair
        (2,3,1,0),  # Both (3) and (2)
        (3,2,0,1),  # Both (4) and (2)
      ])

    for perm in allowed_permutations:
      jirep_perm = np.transpose(jirep, perm)
      maxabs     = np.max(np.abs(jirep_perm - jirep))
      self.assertLessEqual(maxabs, num_zero,
                           msg="Check agreement of eri with permutation " + str(perm)
                              + " failed with difference " + str(maxabs))

  def compare_hf_7_eri_e2e(self, case, hfres):
    num_tol = case["testing"]["numeric_tolerance"]

    if not "eri_ffff" in hfres:
      print("Skipping eri_e2e test for " + case["testing"]["name"]
            + ", since eri_ffff not available")
      return

    noa = hfres["n_orbs_alpha"]
    na  = hfres["n_alpha"]
    nb  = hfres["n_beta"]
    j   = molsturm.compute_coulomb_ff(hfres)
    k   = molsturm.compute_exchange_ff(hfres)

    ene_j = np.sum(j.diagonal()[:na]) + np.sum(j.diagonal()[noa:noa+nb])
    ene_j = 0.5*ene_j
    self.assertAlmostEqual(hfres["energy_coulomb"], ene_j, tol=num_tol,
                           prefix="energy_coulomb from eri: ")

    ene_k = np.sum(k.diagonal()[:na]) + np.sum(k.diagonal()[noa:noa+nb])
    ene_k = -0.5*ene_k
    self.assertAlmostEqual(hfres["energy_exchange"], ene_k, tol=num_tol,
                           prefix="energy_exchange from eri: ")

    e2e = ene_j + ene_k
    hfene = molsturm.compute_derived_hartree_fock_energies(hfres)
    self.assertAlmostEqual(hfene["energy_2e"], e2e, tol=num_tol,
                           prefix="energy_1e from eri: ")


  def __compare_fock_to_orben(self, fock, orben, tolerance, what):
    """Compare a fock matrix to the orben_f values on the
       diagonal

       tolerance:  for comparison
       what:       description what is compared with orben_f
    """
    if fock.shape[0] > 15:
      maxabs = np.max(np.max(fock - np.diag(orben)))
      self.assertLessEqual(maxabs, tolerance,
                           msg="Too large maxabs deviation when comparing "
                           + what + " and orben_f.")
    else:
      self.assertArrayAlmostEqual(fock, np.diag(orben), tol=tolerance,
                                  prefix="Compare " + what + " and orben_f.")


  def compare_hf_8_mos_diagonal(self, case, hfres):
    for k in [ "hcore_ff", "eri_ffff"]:
      if not k in hfres:
        print("Skipping 'fock matrix in mos is diagonal' test for "
              + case["testing"]["name"]
              + ", since " + k + " not available")
        return

    num_tol = case["testing"]["numeric_tolerance"]
    hcore   = hfres["hcore_ff"]
    j       = molsturm.compute_coulomb_ff(hfres)
    k       = molsturm.compute_exchange_ff(hfres)

    fock = hcore + j - k
    self.__compare_fock_to_orben(fock, hfres["orben_f"], tolerance=10*num_tol,
                                 what="fock matrix in MOs")


  def compare_hf_9_fock_diagonal(self, case, hfres):
    if not "fock_ff" in hfres:
      return

    num_tol = case["testing"]["numeric_tolerance"]
    self.__compare_fock_to_orben(hfres["fock_ff"], hfres["orben_f"], tolerance=10*num_tol,
                                 what="fock_ff")


#
# ----------------------------------------------------------------
#

  def compare_hf_results_small(self, case, hfres):
    """
    Compare the most important properties one after another
    """
    self.compare_hf_0_energies(case, hfres)
    self.compare_hf_1_orben_orbcoeff(case, hfres)
    self.compare_hf_2_convergence(case, hfres)
    self.compare_hf_3_spin(case, hfres)


  def compare_hf_results(self, case, hfres):
    """
    Call all compare functions after another
    """
    start_times=dict()
    end_times=dict()
    start_all = time.time()

    for fun in dir(HartreeFockTestCase):
      if fun == "compare_hf_results":
        continue # Do not call ourselves.

      if fun.startswith("compare_hf_"):
        start_times[fun] = time.time()
        getattr(self, fun) (case, hfres)
        end_times[fun] = time.time()

    time_fmt = "   {0:27s}  {1:13.8f}s"
    end_all = time.time()
    if test_constants.PRINT_TIMINGS:
      for k in sorted(start_times):
        print(time_fmt.format(k, end_times[k] - start_times[k]))
      print("   " + 27*"-" + "  " + 13*"-")
      print(time_fmt.format("TOTAL", end_all - start_all)+"\n")

