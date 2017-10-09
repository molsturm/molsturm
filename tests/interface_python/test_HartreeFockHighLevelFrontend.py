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

import molsturm
import unittest
import testdata
from HartreeFockTestCase import HartreeFockTestCase


class TestHartreeFockHighLevelFrontend(HartreeFockTestCase):
    @unittest.skipUnless("sturmian/atomic/cs_static14" in molsturm.available_basis_types,
                         "Required basis type sturmian/atomic/cs_static14 is not "
                         "available")
    def test_sturmian_basis_on_the_fly(self):
        case = testdata.test_cases_by_name("be_cs32")[0]

        try:
            params = molsturm.ScfParameters.from_dict(case["input_parameters"])
        except (ValueError, KeyError, TypeError) as e:
            raise unittest.SkipTest("Skipped subtest " + case["testing"]["name"] +
                                    ", since construction of ScfParameters failed: " +
                                    str(e))
        system = params.system
        discretisation = params["discretisation"]
        scf = params["scf"]
        hfres = molsturm.hartree_fock(system, basis_type=discretisation["basis_type"],
                                      conv_tol=scf["max_error_norm"],
                                      max_iter=scf["max_iter"],
                                      n_eigenpairs=scf["n_eigenpairs"],
                                      restricted=(scf["kind"] == "restricted"),
                                      guess=params["guess/method"],
                                      eigensolver=scf["eigensolver/method"],
                                      print_iterations=scf["print_iterations"],
                                      k_exp=discretisation["k_exp"],
                                      n_max=discretisation["n_max"],
                                      l_max=discretisation["l_max"],
                                      m_max=discretisation["m_max"])

        # TODO test guess from previous via this interface, too

        self.compare_hf_results_small(case, hfres)

    @unittest.skipUnless("gaussian/libint" in molsturm.available_basis_types,
                         "Required basis type gaussian/libint is not available")
    def test_gaussian_basis_on_the_flay(self):
        case = testdata.test_cases_by_name("c_321g")[0]

        try:
            params = molsturm.ScfParameters.from_dict(case["input_parameters"])
        except (ValueError, KeyError, TypeError) as e:
            raise unittest.SkipTest("Skipped subtest " + case["testing"]["name"] +
                                    ", since construction of ScfParameters failed: " +
                                    str(e))
        system = params.system
        discretisation = params["discretisation"]
        scf = params["scf"]
        hfres = molsturm.hartree_fock(system, basis_type=discretisation["basis_type"],
                                      conv_tol=scf["max_error_norm"],
                                      max_iter=scf["max_iter"],
                                      n_eigenpairs=scf["n_eigenpairs"],
                                      restricted=(scf["kind"] == "restricted"),
                                      guess=params["guess/method"],
                                      eigensolver=scf["eigensolver/method"],
                                      print_iterations=scf["print_iterations"],
                                      basis_set_name=discretisation["basis_set_name"])

        # TODO test guess from previous via this interface, too

        self.compare_hf_results_small(case, hfres)
