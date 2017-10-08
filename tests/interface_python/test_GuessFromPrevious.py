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
from molsturm.scf_guess import extrapolate_from_previous
import unittest
import testdata
from HartreeFockTestCase import HartreeFockTestCase

class TestFromPrevious(HartreeFockTestCase):
  def __run_test(self, name):
    case = testdata.test_cases_by_name(name)[0]

    try:
      inp = case["input_parameters"]
      scfparams = molsturm.ScfParameters.from_dict(inp)

      inp.setdefault("scf", {})
      inp["scf"]["conv_tol"] = 1e-5
      scfparams_lowtol = molsturm.ScfParameters.from_dict(inp)
    except (ValueError, KeyError, TypeError) as e:
      raise unittest.SkipTest("Skipped subtest " + testing["name"] + ", since"
                              "construction of ScfParameters failed: " + str(e))

    res = molsturm.self_consistent_field(scfparams_lowtol)

    # Construct extrapolated guess and run again
    guess = extrapolate_from_previous(res, scfparams)
    scfparams.set_guess_external(*guess)
    hfres = molsturm.self_consistent_field(scfparams)

    self.assertLessEqual(hfres["final_error_norm"], res["final_error_norm"])
    self.assertLessEqual(hfres["n_iter"], 3)
    self.compare_hf_results_small(case, hfres)


  def test_restricted_be_cs32(self):
    self.__run_test("be_cs32")

  def test_unrestricted_c_cs41(self):
    self.__run_test("c_cs41")

  def test_unrestricted_c_sto3g(self):
    self.__run_test("c_321g")
