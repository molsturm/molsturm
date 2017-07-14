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

from FciTestCase import FciTestCase
import testdata
import molsturm.posthf
from molsturm import INPUT_PARAMETER_KEY

class TestFCI(FciTestCase):
  """This test should assure that molsturm results stay the same
     between code changes or algorithm updates
  """

  @classmethod
  def setUpClass(cls):
    cls.cases = testdata.test_cases()

  def test_fci(self):
    for case in self.cases:
      testing   = case["testing"]
      hf        = case["hf"]

      try:
        fciparams = case["fci"][INPUT_PARAMETER_KEY]
      except KeyError:
        # No fci test data available for this case
        continue

      fci = molsturm.posthf.fci(hf, **fciparams)
      self.compare_fci_results(case, fci)

