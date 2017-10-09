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


class FciTestCase(NumCompTestCase):
    def compare_fci_0_energies(self, case, fci):
        num_tol = case["testing"]["numeric_tolerance"]
        ref_states = case["fci"]["states"]
        fci_states = fci["states"]

        self.assertEqual(len(fci_states), len(ref_states))
        for i in range(len(fci_states)):
            self.assertAlmostEqual(fci_states[i]["energy"], ref_states[i]["energy"],
                                   tol=num_tol,
                                   prefix="root " + str(i) + " energy: ")

    def compare_fci_1_spin(self, case, fci):
        num_tol = case["testing"]["numeric_tolerance"]
        ref_states = case["fci"]["states"]
        fci_states = fci["states"]

        self.assertEqual(len(fci_states), len(ref_states))
        for i in range(len(fci_states)):
            for key in ["multiplicity", "spin_squared"]:
                if ref_states[i][key] is None:
                    print("Skipping FCI test for " + key + " for " +
                          case["testing"]["name"] +
                          ", since no reference data available.")
                    continue  # Skip keys where the reference value is not available.

                self.assertAlmostEqual(fci_states[i][key], ref_states[i][key],
                                       tol=num_tol,
                                       prefix="root " + str(i) + " " + key + " ")

    def compare_fci_results(self, case, fci):
        """
        Call all compare functions after another
        """
        for fun in dir(FciTestCase):
            if fun == "compare_fci_results":
                continue  # Do not call ourselves.

            if fun.startswith("compare_fci_"):
                getattr(self, fun)(case, fci)
