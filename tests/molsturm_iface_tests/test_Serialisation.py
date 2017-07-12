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

import molsturm
import numpy as np
import os
import tempfile
import testdata
import unittest

class TestSerialisation(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.cases = testdata.test_cases()

  def assert_equal(self,d1,d2):
    self.assertListEqual([d1.keys()],[d2.keys()])

    for k in d1:
      if isinstance(d1[k],np.ndarray):
        self.assertTrue( np.all(d1[k] == d2[k]) )
      elif isinstance(d1[k],dict):
        self.assert_equal(d1[k],d2[k])
      elif isinstance(d1[k],list):
        self.assertListEqual(d1[k],d2[k])
      else:
        self.assertEqual(type(d1[k]), type(d2[k]))
        self.assertEqual(d1[k],d2[k])


  # ----------------------------------------------------------

  def test_yaml_serialisation(self):
    for case in self.cases:
      with self.subTest(label=case["testing"]["name"]):
        hfres = case["hf"]

        tmp = tempfile.mktemp(suffix=".yaml")
        molsturm.dump_yaml(hfres, tmp)
        back = molsturm.load_yaml(tmp)
        os.remove(tmp)

        self.assert_equal(hfres, back)

  def test_hdf5_serialisation(self):
    for case in self.cases:
      with self.subTest(label=case["testing"]["name"]):
        hfres = case["hf"]

        tmp = tempfile.mktemp(suffix=".hdf5")
        molsturm.dump_hdf5(hfres, tmp)
        back = molsturm.load_hdf5(tmp)
        os.remove(tmp)

        self.assert_equal(hfres, back)

