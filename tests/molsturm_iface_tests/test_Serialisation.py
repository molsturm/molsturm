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
import molsturm
import tempfile
import os
import numpy as np

class TestSerialisation(unittest.TestCase):
  @classmethod
  def import_cases(cls):
    dir_of_this_script = os.path.dirname( os.path.abspath( __file__ ) )

    cls.cases = dict()

    for f in os.listdir(dir_of_this_script+"/testdata"):
      fn,ext = os.path.splitext(f)

      # Ignore input files and only consider yaml files
      if ext != ".yaml": continue
      if fn[-3:] == ".in": continue

      key = os.path.basename(fn)
      datafile = dir_of_this_script+"/testdata/"+key+".yaml"
      cls.cases[key] = molsturm.load_yaml(datafile)

  @classmethod
  def setUpClass(cls):
    cls.import_cases()

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
        self.assertEqual(d1[k],d2[k])


  # ----------------------------------------------------------

  def test_yaml_serialisation(self):
    for caselal in self.cases:
      with self.subTest(label=caselal):
        params = self.cases[caselal]

        tmp = tempfile.mktemp(suffix=".yaml")
        molsturm.dump_yaml(params,tmp)
        back = molsturm.load_yaml(tmp)
        os.remove(tmp)

        self.assert_equal(params,back)

  def test_hdf5_serialisation(self):
    for caselal in self.cases:
      with self.subTest(label=caselal):
        params = self.cases[caselal]

        tmp = tempfile.mktemp(suffix=".yaml")
        molsturm.dump_hdf5(params,tmp)
        back = molsturm.load_hdf5(tmp)
        os.remove(tmp)

        self.assert_equal(params,back)

