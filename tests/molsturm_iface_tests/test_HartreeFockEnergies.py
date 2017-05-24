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
import data_hf_energies as data
import numpy as np

class TestHartreeFockEnergies(NumCompTestCase):
  """This test should ensure, that we our molsturm can reproduce
     data exactly in the way computed with ORCA or another standard
     quantum chemistry program.
  """

  def run_case(self, case):
    inparams = dict(data.params_all)
    inparams.update(case["system"])
    inparams.update(case["params"])

    hf = molsturm.hartree_fock(**inparams)
    hf_1e = hf["energy_kinetic"] + hf["energy_nuclear_attraction"]
    hf_2e = hf["energy_coulomb"] + hf["energy_exchange"]

    error = inparams["error"]
    ref = case["reference"]

    self.assertAlmostEqual(hf_1e, ref["energy_1e"], tol=10*error, prefix="1e energy: ")
    self.assertAlmostEqual(hf_2e, ref["energy_2e"], tol=10*error, prefix="2e energy: ")
    self.assertAlmostEqual(hf["energy_nuclear_repulsion"], ref["energy_nucrep"], 
                           tol=10*error, prefix="Nuclear repulsion: ")
    self.assertAlmostEqual(hf["energy_total"], ref["energy_total"], tol=error,
                           prefix="Total energy: ")

    if "energy_mos" in ref:
      tol=max(1e-6, error)   # Unfortunately ORCA spits out the mos only up to very low accuracy
      mos = np.array(ref["energy_mos"])
      self.assertArrayAlmostEqual(hf["orbital_energies_f"], mos, tol=tol, prefix="MO energies: ")

  def test_run_all_cases(self):
    for case in data.cases:
      with self.subTest(msg=case["description"]):
        self.run_case(case)

