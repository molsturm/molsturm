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
import unittest
import testdata
from HartreeFockTestCase import HartreeFockTestCase

class TestFromPrevious(HartreeFockTestCase):
  def __run_test(self, name):
    case = testdata.test_cases_by_name(name)[0]
    params = case["params"]

    if not params["basis_type"] in molsturm.available_basis_types:
      raise unittest.SkipTest("Required basis type "
                              + params["basis_type"] + " is not available")

    # Remove tolerance
    params_notol = dict(params)
    if "conv_tol" in params_notol:
      del params_notol["conv_tol"]

    mintol = 1e-5
    res = molsturm.hartree_fock(**params_notol, conv_tol=mintol)

    if "guess" in params:
      params_copy=dict(params)
      del params_copy["guess"]
      params = params_copy

    hfres = molsturm.hartree_fock(**params, guess=res )

    self.assertLessEqual(hfres["final_error_norm"], res["final_error_norm"])
    self.assertLessEqual(hfres["n_iter"], 3)
    self.compare_hf_results_small(case, hfres)


  def test_restricted_be_cs32(self):
    self.__run_test("be_cs32")

  def test_unrestricted_c_cs41(self):
    self.__run_test("c_cs41")

  def test_unrestricted_c_sto3g(self):
    self.__run_test("c_sto3g")
