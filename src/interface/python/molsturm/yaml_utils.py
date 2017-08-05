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

import numpy as np
import yaml

def ndarray_representer(dumper, data):
  """YAML representer for numpy arrays"""
  return dumper.represent_sequence('!ndarray', data.tolist())


def numpy_scalar_representer(dumper, data):
  """YAML representer for numpy scalar values.
  They are represented as their python equivalents.
  """
  return dumper.represent_data(np.asscalar(data))


def ndarray_constructor(loader, node):
  """YAML constructor for numpy arrays"""
  value = loader.construct_sequence(node)
  return np.array(value)


def install_constructors():
  """Install all YAML constructors defined in this module"""
  yaml.constructor.SafeConstructor.add_constructor("!ndarray", ndarray_constructor)


def install_representers():
  """Install all YAML representers defined in this module"""
  yaml.representer.SafeRepresenter.add_representer(np.ndarray, ndarray_representer)

  for np_type_category in [ 'complex', 'float', 'int' ]:
    for tpe in np.sctypes[np_type_category]:
      yaml.representer.SafeRepresenter.add_representer(tpe, numpy_scalar_representer)
