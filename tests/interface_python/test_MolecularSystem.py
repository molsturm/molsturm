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
from molsturm.MolecularSystem import MolecularSystem

@unittest.skipUnless("sturmian/atomic/cs_static14" in molsturm.available_basis_types,
                     "Required basis type sturmian/atomic/cs_static14 is not available")
class TestFromPrevious(HartreeFockTestCase):
  def test_dummy(self):
    case = testdata.test_cases_by_name("be_cs32")[0]
    params = case["params"]
    params_copy = dict(params)

    system = MolecularSystem(**{ k : params.get(k,None) for k in MolecularSystem._fields })
    for k in MolecularSystem._fields:
      if k in params_copy:
        del params_copy[k]

    hfres = molsturm.hartree_fock(system, **params_copy)
    self.compare_hf_results_small(case, hfres)

