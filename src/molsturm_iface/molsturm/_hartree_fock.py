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

import molsturm_iface as iface
import numpy as np
from ._constants import HFRES_ARRAY_KEYS, HFRES_INPUT_PARAMETER_KEY
from collections import Iterable

def __to_double_vector(val):
  ret = iface.DoubleVector()
  for v in val.flatten():
    ret.push_back(float(v))
  return ret


def __to_coords(arr):
  if not isinstance(arr,Iterable):
    raise TypeError("Argument passed to \"coords\" needs to be iterable")

  ret = iface.DoubleVector()
  def parse_coord(c):
    if len(c) != 3:
      raise ValueError("All items of the coordinates array need have exactly 3 items.")
    for i in range(3):
      try:
        ret.push_back(float(c[i]))
      except ValueError:
        raise ValueError("All elements of coord vectors need to be castible to float")

  if isinstance(arr[0], Iterable):
    for c in arr:
      parse_coord(c)
  else:
    parse_coord(arr)
  return ret


def __to_nlm(arr):
  ret = iface.IntVector()
  for c in arr:
    if len(c) != 3:
      raise ValueError("Nlm array needs to be of shape n times 3")
    for i in range(3):
      ret.push_back(c[i])
  return ret


def __to_atom_numbers(li):
  ret = iface.IntVector()
  if isinstance(li,int):
    ret.push_back(int(li))
  elif isinstance(li,Iterable):
    for n in li:
      ret.push_back(int(n))
  else:
    raise ValueError("atom numbers needs be an int or a list of ints")
  return ret


def __to_atoms(li):
  ret = iface.StringVector()
  if isinstance(li,str):
    ret.push_back(str(li))
  elif isinstance(li,Iterable):
    for n in li:
      ret.push_back(str(n))
  else:
    raise ValueError("atoms needs to be a string or a list of strings")
  return ret

# Special input parameters for which the above conversion functions need
# to be used before assignment
__params_transform_maps = {
  "coords":       __to_coords,
  "atom_numbers": __to_atom_numbers,
  "atoms":        __to_atoms,
  "nlm_basis":    __to_nlm,
}

# These keys exist in the Parameters struct, but should not be exposed
# via the hartree_fock interface
HF_PARAMS_EXCLUDE_KEYS = [ "restricted_set_by_user", "all" ]

"""The list of keys understood by the hartree_fock function"""
hartree_fock_keys = [ k for k in dir(iface.Parameters)
                      if k[0] != "_" and not k in HF_PARAMS_EXCLUDE_KEYS ]

def hartree_fock(forward_parameters=True, **kwargs):
  """
  Run a Hartree-Fock calculation with molsturm. The list of valid input
  parameters can be retrieved by the means of the list "hartree_fock_keys".

  Parameters:
    forward_parameters:   Should the list of kwargs be forwarded to the
                          output dictionary for archiving purposes?
                          If True(default) the returned dictionary
                          will contain a key "input_parameters" which
                          is a set of parameters to reproduce the very
                          computation. Some transformations may have been
                          applied to the input kwargs (e.g. xyz files
                          are read and stored as a list of atoms and
                          coordinates)

  A couple of selected kwargs:
    - basis_type        The type of the basis used for the calculation
    - coords            List of iterables of size 3: Coordinats of the atoms
    - atoms             List of the atom symbols (in the same order as coords)
  """
  # TODO Better return a dict-like class instead of a dict. That way
  #      we can use the class more easily and distinguish between results
  #      at different levels better.

  # The list of valid keys is the list of keys
  # with the special ones (starting with __) removed.
  params_keys = [ k for k in dir(iface.Parameters) if k[0] != "_" ]
  res_keys = [ k for k in dir(iface.HfResults) if k[0] != "_" ]

  # Build params and run:
  params = iface.Parameters()
  for key in kwargs:
    if not key in params_keys or key in HF_PARAMS_EXCLUDE_KEYS:
      raise ValueError("Keyword " + key + " is unknown to hartree_fock")
    elif key in __params_transform_maps:
      setattr(params,key,__params_transform_maps[key](kwargs[key]))
    else:
      setattr(params,key,kwargs[key])

  # Make a note that the user specified the restricted keyword
  if "restricted" in kwargs: params.restricted_set_by_user = True

  res = iface.hartree_fock(params)

  # Build output dictionary:
  shape_lookup = { "f": res.n_orbs_alpha + res.n_orbs_beta,
                   "b": res.n_bas }

  out = { k :getattr(res,k) for k in res_keys if not k in HFRES_ARRAY_KEYS }
  for k in HFRES_ARRAY_KEYS:
    # Build the shape to cast the numpy arrays into from the
    # suffixes (e.g. _ffff, _bf) and the shape lookup object
    # we created above
    target_shape = tuple( shape_lookup[c] for c in k[k.rfind("_")+1:] )
    ary = np.array(getattr(res,k))
    if ary.size != 0:
      # If the size is 0, then the data has not been computed,
      # so we can ignore it
      out[k] = ary.reshape(target_shape)

  if forward_parameters:
    out[HFRES_INPUT_PARAMETER_KEY] = kwargs

  return out

