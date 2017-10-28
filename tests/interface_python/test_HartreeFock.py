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

from HartreeFockTestCase import HartreeFockTestCase
from testdata import predicates
import molsturm
import test_constants
import testdata
import unittest


class TestHartreeFock(HartreeFockTestCase):
    """This test should assure that molsturm results stay the same
       between code changes or algorithm updates
    """
    @classmethod
    def setUpClass(cls):
        # Select predicate to filter test cases:
        if test_constants.RUN_EXPENSIVE:
            pred = predicates.all()
        else:
            pred = predicates.is_not_expensive()
        cls.cases = testdata.test_cases_by_pred(pred)

    def test_hf(self):
        for case in self.cases:
            testing = case["testing"]

            with self.subTest(label=testing["name"]):
                try:
                    scfparams = molsturm.ScfParameters.from_dict(case["input_parameters"])
                    scfparams.normalise()
                except (ValueError, KeyError, TypeError) as e:
                    raise unittest.SkipTest("Skipped subtest " + testing["name"] +
                                            ", since construction of ScfParameters "
                                            "failed: " + str(e))

                hfres = molsturm.self_consistent_field(scfparams)
                self.compare_hf_results(case, hfres)
