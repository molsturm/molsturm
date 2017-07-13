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

from NumCompTestCase import NumCompTestCase

class FciTestCase(NumCompTestCase):
  def compare_fci_0_energies(self, case, fci):
    num_tol = case["testing"]["numeric_tolerance"]
    ref     = case["fci"]

    for i in range(len(fci)):
      self.assertAlmostEqual(fci[i]["energy"], ref[i]["energy"], tol=num_tol,
                             prefix="root " + str(i) + " energy: ")


  def compare_fci_results(self, case, fci):
    """
    Call all compare functions after another
    """
    for fun in dir(FciTestCase):
      if fun == "compare_fci_results":
        continue # Do not call ourselves.

      if fun.startswith("compare_fci_"):
        getattr(self, fun) (case, fci)



