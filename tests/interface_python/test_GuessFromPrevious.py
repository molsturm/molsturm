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
from molsturm.scf_guess import extrapolate_from_previous
import unittest
import testdata
from HartreeFockTestCase import HartreeFockTestCase


class TestFromPrevious(HartreeFockTestCase):
    def __run_test_same(self, name):
        case = testdata.test_cases_by_name(name)[0]

        inp = case["input_parameters"]
        scfparams = molsturm.ScfParameters.from_dict(inp)
        scfparams_lowtol = molsturm.ScfParameters.from_dict(inp)
        scfparams_lowtol["scf/con_tol"] = 1e-5

        try:
            scfparams.normalise()
            scfparams_lowtol.normalise()
        except (ValueError, KeyError, TypeError) as e:
            raise unittest.SkipTest("Skipped subtest " + case["testing"]["name"] +
                                    ", since construction of ScfParameters "
                                    "failed: " + str(e))

        res = molsturm.self_consistent_field(scfparams_lowtol)

        # Construct extrapolated guess and run again
        guess = extrapolate_from_previous(res, scfparams)
        scfparams.set_guess_external(*guess)
        hfres = molsturm.self_consistent_field(scfparams)

        self.assertLessEqual(hfres["final_error_norm"], res["final_error_norm"])
        self.assertLessEqual(hfres["n_iter"], 3)
        self.compare_hf_results_small(case, hfres)

    def __run_test_sturmian_interpolate(self, name, n_max, l_max, m_max):
        case = testdata.test_cases_by_name(name)[0]

        inp = case["input_parameters"]
        scfparams = molsturm.ScfParameters.from_dict(inp)

        # This done to stabilise the outcome such that we don't
        # get the correct minimum once and a wrong minimum at another time
        scfparams["guess/method"] = "random"

        try:
            scfparams.normalise()
        except (ValueError, KeyError, TypeError) as e:
            raise unittest.SkipTest("Skipped subtest " + case["testing"]["name"] +
                                    ", since construction of ScfParameters "
                                    "failed: " + str(e))

        smallparams = scfparams.copy()
        origbasis = scfparams.basis

        smallparams.basis = molsturm.Basis.construct(
            "sturmian/atomic/" + origbasis.backend,
            scfparams.system, n_max=n_max, l_max=l_max,
            m_max=m_max, k_exp=origbasis.k_exp,
        )
        smallres = molsturm.self_consistent_field(smallparams)

        res_no_guess = molsturm.self_consistent_field(scfparams)

        guess = molsturm.scf_guess.extrapolate_from_previous(smallres, scfparams)
        scfparams.set_guess_external(*guess)
        res_guess = molsturm.self_consistent_field(scfparams)

        self.assertAlmostEqual(res_no_guess["energy_ground_state"],
                               res_guess["energy_ground_state"],
                               tol=case["testing"]["numeric_tolerance"],
                               prefix="Energy of ground state without and with guess: ")

        self.assertLess(res_guess["n_iter"], res_no_guess["n_iter"],
                        msg="Number of iterations with guess not smaller than without.")

    def test_restricted_be_cs32(self):
        self.__run_test_same("be_cs32")
        self.__run_test_sturmian_interpolate("be_cs32", n_max=3, l_max=1, m_max=1)

    def test_unrestricted_c_cs41(self):
        self.__run_test_same("c_cs41")
        self.__run_test_sturmian_interpolate("c_cs41", n_max=3, l_max=1, m_max=1)

    def test_unrestricted_c_sto3g(self):
        self.__run_test_same("c_321g")
