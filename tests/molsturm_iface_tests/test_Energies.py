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
import molsturm.posthf
import numpy as np
import unittest

import data_energies as data

class TestEnergies(NumCompTestCase):
  """This test should ensure, that we our molsturm can reproduce
     data exactly in the way computed with ORCA or another standard
     quantum chemistry program.
  """

  @classmethod
  def import_cases(cls):
    # Parse input and build test cases:
    return data.cases

  @classmethod
  def setUpClass(cls):
    cls.cases = cls.import_cases()
    cls.hf_results = dict()

  def run_hf(self,case):
    error = case["params"]["error"]

    do_postHf = False
    if "mp2" in molsturm.posthf.available_methods:
      do_postHf = True

    params = dict(case["params"])
    params.update({
      "eigensolver":    "lapack",
      "guess_esolver":  "lapack",
      #
      "export_repulsion_integrals": do_postHf,
    })

    hf = molsturm.hartree_fock(**params)
    hf_1e = hf["energy_kinetic"] + hf["energy_nuclear_attraction"]
    hf_2e = hf["energy_coulomb"] + hf["energy_exchange"]
    ref = case["reference_hf"]

    self.assertAlmostEqual(hf_1e, ref["energy_1e"], tol=10*error, prefix="1e energy: ")
    self.assertAlmostEqual(hf_2e, ref["energy_2e"], tol=10*error, prefix="2e energy: ")
    self.assertAlmostEqual(hf["energy_nuclear_repulsion"], ref["energy_nucrep"], 
                           tol=10*error, prefix="Nuclear repulsion: ")
    self.assertAlmostEqual(hf["energy_total"], ref["energy_total"], tol=error,
                           prefix="Total energy: ")

    if "energy_mos" in ref:
      tol=max(1e-6, error)   # Unfortunately ORCA spits out the mos only up to very low accuracy
      mos = np.array(ref["energy_mos"])
      self.assertArrayAlmostEqual(hf["orben_f"], mos, tol=tol, prefix="MO energies: ")
    return hf

  def run_mp2(self,case,hf):
    if not "reference_mp2" in case: return None

    error = case["params"]["error"]

    mp2 = molsturm.posthf.mp2(hf, max_memory=64*1024*1024)
    mp2_ref = case["reference_mp2"]
    self.assertAlmostEqual(mp2["energy_mp2"], mp2_ref["energy_mp2"],
                           tol=error, prefix="MP2 energy: ")
    self.assertAlmostEqual(mp2["energy_ground_state"], mp2_ref["energy_ground_state"],
                           tol=error, prefix="Ground state energy: ")

  # --------------------------------------------

  def test_0_hf(self):
    for caselal in self.cases:
      with self.subTest(label=caselal):
        self.hf_results[caselal] = self.run_hf(self.cases[caselal])

  @unittest.skipUnless("mp2" in molsturm.posthf.available_methods,
                      "mp2 not available => Skipping mp2 tests")
  def test_1_mp2(self):
    for caselal in self.cases:
      with self.subTest(label=caselal):
        if not caselal in self.hf_results:
          raise self.fail("HF results not available")
        else:
          self.run_mp2(self.cases[caselal],self.hf_results[caselal])

