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

from molsturm import INPUT_PARAMETER_KEY
from MP2TestCase import MP2TestCase
from testdata import predicates
import molsturm.posthf
import testdata
import unittest

# TODO Note that this class is very similar to test_FCI and could
#      perhaps be unified with it


@unittest.skipUnless("mp2" in molsturm.posthf.available_methods,
                     "mp2 not available => Skipping mp2 tests")
class TestMP2(MP2TestCase):
    """This test should assure that molsturm results stay the same
       between code changes or algorithm updates
    """

    @classmethod
    def setUpClass(cls):
        cls.cases = testdata.test_cases_by_pred(predicates.tests_method("mp2"))

    def test_mp2(self):
        for case in self.cases:
            hf = case["hf"]
            mp2params = case["mp2"][INPUT_PARAMETER_KEY]

            mp2 = molsturm.posthf.mp2(hf, **mp2params)
            self.compare_mp2_results(case, mp2)
