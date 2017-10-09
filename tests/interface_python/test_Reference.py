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
                scfparams = molsturm.ScfParameters.from_dict(case["input_parameters"])

                # Update parameters if posthf is done
                hf = molsturm.self_consistent_field(scfparams)

                self.compare_hf_results(case, hf)
                self.hf_results[testing["name"]] = hf

    @unittest.skipUnless("mp2" in molsturm.posthf.available_methods,
                         "mp2 not available => Skipping mp2 tests")
    def test_1_mp2(self):
        for case in self.cases:
            if "mp2" not in case:
                continue

            testing = case["testing"]
            name = testing["name"]
            params = case["mp2"][INPUT_PARAMETER_KEY]
            with self.subTest(label=name):
                if name not in self.hf_results:
                    raise self.fail("HF results not available")

                mp2 = molsturm.posthf.mp2(self.hf_results[name], **params)
                self.compare_mp2_results(case, mp2)

    @unittest.skipUnless("fci" in molsturm.posthf.available_methods,
                         "fci not available => Skipping fci tests")
    def test_2_fci(self):
        for case in self.cases:
            if "fci" not in case:
                continue

            testing = case["testing"]
            name = testing["name"]
            params = case["fci"][INPUT_PARAMETER_KEY]
            with self.subTest(label=name):
                if name not in self.hf_results:
                    raise self.fail("HF results not available")

                fci = molsturm.posthf.fci(self.hf_results[name], **params)
                self.compare_fci_results(case, fci)
