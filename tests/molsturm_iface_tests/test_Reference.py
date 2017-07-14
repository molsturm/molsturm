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

from FciTestCase import FciTestCase
from HartreeFockTestCase import HartreeFockTestCase
from MP2TestCase import MP2TestCase
from molsturm import INPUT_PARAMETER_KEY
import molsturm
import molsturm.posthf
import testdata
import unittest

@unittest.skipUnless("gaussian/libint" in molsturm.available_basis_types,
                     "gaussian/libint not available. Skipping reference tests")
class TestReference(HartreeFockTestCase, MP2TestCase, FciTestCase):
  """This test should ensure, that we our molsturm can reproduce
     data exactly in the way computed with ORCA or another standard
     quantum chemistry program.
  """

  @classmethod
  def setUpClass(cls):
    cls.cases = testdata.reference_cases()
    cls.hf_results = dict()

  # --------------------------------------------

  def test_0_hf(self):
    for case in self.cases:

      testing = case["testing"]
      name = testing["name"]
      with self.subTest(label=name):
        params = case["params"]

        # Update parameters if posthf is done
        if testing["any_posthf"]:
          params.setdefault("export_fock_matrix",         True)
          params.setdefault("export_repulsion_integrals", True)
          params.setdefault("export_hcore_matrix",        True)
        hf = molsturm.hartree_fock(**params)

        self.compare_hf_results(case, hf)
        self.hf_results[ testing["name"] ] = hf


  @unittest.skipUnless("mp2" in molsturm.posthf.available_methods,
                       "mp2 not available => Skipping mp2 tests")
  def test_1_mp2(self):
    for case in self.cases:
      if not "mp2" in case:
        continue

      testing = case["testing"]
      name = testing["name"]
      params = case["mp2"][INPUT_PARAMETER_KEY]
      with self.subTest(label=name):
        if not name in self.hf_results:
          raise self.fail("HF results not available")

        mp2 = molsturm.posthf.mp2(self.hf_results[name], **params)
        self.compare_mp2_results(case, mp2)


  @unittest.skipUnless("fci" in molsturm.posthf.available_methods,
                       "fci not available => Skipping fci tests")
  def test_2_fci(self):
    for case in self.cases:
      if not "fci" in case:
        continue

      testing = case["testing"]
      name = testing["name"]
      params = case["fci"][INPUT_PARAMETER_KEY]
      with self.subTest(label=name):
        if not name in self.hf_results:
          raise self.fail("HF results not available")

        fci = molsturm.posthf.fci(self.hf_results[name], **params)
        self.compare_fci_results(case, fci)

