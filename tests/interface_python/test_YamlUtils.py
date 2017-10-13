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

import unittest
import numpy as np
import yaml
import molsturm


class TestYamlUtils(unittest.TestCase):
    def __test_yaml_array(self, array):
        molsturm.yaml_utils.install_representers()
        molsturm.yaml_utils.install_constructors()
        s = yaml.safe_dump(array)
        res = yaml.safe_load(s)

        self.assertTrue(np.all(array == res), msg="Original array (" +
                        str(array) + ") and represented array(" +
                        str(res) + ") not equivalent.")

    def test_yaml_conrep_1darray(self):
        self.__test_yaml_array(np.array([1, 2, 3, 4]))

    def test_yaml_conrep_2darray(self):
        arr = [[1, 2], [3, 4], [5, 6]]
        arr = np.array(arr)
        self.__test_yaml_array(arr)

    def test_yaml_conrep_3darray(self):
        arr = [[[1, 2], [3, 4], [5, 6]],
               [[7, 8], [9, 0], [1, 2]]]
        arr = np.array(arr)
        self.__test_yaml_array(arr)
