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

from NumCompTestCase import NumCompTestCase
import collections
import molsturm
import molsturm.sturmian
import unittest


class TestSturmian(NumCompTestCase):
    def __run_find_kopt(self, case):
        scfparams = molsturm.ScfParameters()
        scfparams.system = molsturm.MolecularSystem(case.atom,
                                                    multiplicity=case.multiplicity)

        if case.n_max == 3 and case.l_max == 2:
            basis_type = "sturmian/atomic/cs_static14"
        else:
            basis_type = "sturmian/atomic/cs_reference_pc"

        if basis_type not in molsturm.available_basis_types:
            raise unittest.SkipTest("Skipped subtest " + case.atom +
                                    ", since basis type not available: " +
                                    basis_type)

        scfparams.basis = molsturm.Basis.construct(basis_type,
                                                   scfparams.system,
                                                   k_exp=0.0,  # dummy
                                                   n_max=case.n_max,
                                                   l_max=case.l_max)

        scfparams["scf/eigensolver/method"] = "lapack"
        scfparams["guess/eigensolver/method"] = "lapack"

        conv_tol = 1e-3
        best = molsturm.sturmian.cs.find_kopt(scfparams, conv_tol=conv_tol)
        k = best["input_parameters"]["discretisation"]["k_exp"]

        self.assertAlmostEqual(k, case.ref_kexp, tol=conv_tol,
                               prefix="optimal exponent kopt: ")
        self.assertAlmostEqual(best["energy_ground_state"], case.ref_energy,
                               tol=conv_tol / 100,
                               prefix="ground state energy: ")

    def test_find_kopt(self):
        KoptTestCase = collections.namedtuple("KoptTestCase",
                                              ["atom", "multiplicity",
                                               "n_max", "l_max",
                                               "ref_kexp", "ref_energy"])

        cases = [KoptTestCase("C",  3, 3, 2, 2.926, -36.21728187346766),
                 KoptTestCase("Be", 1, 3, 2, 2.032, -14.118626126201944),
                 KoptTestCase("Ne", 1, 5, 1, 4.585, -128.0943438001784)]

        for case in cases:
            with self.subTest(label=case.atom):
                self.__run_find_kopt(case)
