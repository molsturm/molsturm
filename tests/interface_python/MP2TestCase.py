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


class MP2TestCase(NumCompTestCase):
    def compare_mp2_0_energies(self, case, mp2):
        num_tol = case["testing"]["numeric_tolerance"]
        ref = case["mp2"]

        for key in ref:
            if not key.startswith("energy_"):
                continue
            self.assertAlmostEqual(mp2[key], ref[key], tol=num_tol, prefix=key + ": ")

    def compare_mp2_results(self, case, mp2):
        """
        Call all compare functions after another
        """
        for fun in dir(MP2TestCase):
            if fun == "compare_mp2_results":
                continue  # Do not call ourselves.

            if fun.startswith("compare_mp2_"):
                getattr(self, fun)(case, mp2)
