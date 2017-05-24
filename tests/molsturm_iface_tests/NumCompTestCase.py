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

import unittest
import numpy as np

class NumCompTestCase(unittest.TestCase):
  def assertAlmostEqual(self,lhs,rhs,tol,prefix=""):
    self.assertTrue(np.isclose(lhs,rhs,rtol=tol,atol=tol),
                    msg=prefix + "lhs ("+str(lhs) + ") not equal to rhs ("+str(rhs)+")"
                    " (tolerance == " + str(tol) + ")")

  def assertArrayAlmostEqual(self,lhs,rhs,tol,prefix=""):
    self.assertEqual(lhs.shape,rhs.shape,
                     msg=prefix+"lhs shape ("+str(lhs.shape)+") different from "
                     "rhs shape ("+str(rhs.shape)+").")

    rhsf = rhs.flatten()
    lhsf = lhs.flatten()
    for i in range(len(rhsf)):
      self.assertTrue(np.isclose(lhsf[i],rhsf[i],rtol=tol,atol=tol),
                      msg=prefix + "element "+str(i) +" of lhs and lhs differ:\n"+\
                      str(lhsf[i]) + " != " + str(rhsf[i]) + " (tolerance: " + str(tol) +\
                      "). lhs array:\n" + str(lhs) + "\nrhs array:\n"+str(rhs))

