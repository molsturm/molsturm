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

def __np_to_coords(arr):
  ret = iface.DoubleVector()
  for c in arr:
    if len(c) != 3:
      raise ValueError("All items of the coordinates array need have exactly 3 items.")
    for i in range(3):
      ret.push_back(c[i])
  return ret

def __np_to_nlm(arr):
  ret = iface.IntVector()
  for c in arr:
    if len(c) != 3:
      raise ValueError("Nlm array needs to be of shape n times 3")
    for i in range(3):
      ret.push_back(c[i])
  return ret

def __to_int_vector(li):
  ret = iface.IntVector()
  for n in li:
    ret.push_back(n)
  return ret

def __to_string_vector(li):
  ret = iface.StringVector()
  for n in li:
    ret.push_back(n)
  return ret

# Special input parameters for which the above conversion functions need
# to be used before assignment
__params_transform_maps = {
  "coords":       __np_to_coords,
  "atom_numbers": __to_int_vector,
  "atoms":        __to_string_vector,
  "nlm_basis":    __np_to_nlm,
}

"""The list of keys understood by the hartree_fock function"""
hartree_fock_keys = [ k for k in dir(iface.Parameters) if k[0] != "_" ]

def hartree_fock(forward_parameters=True, **kwargs):
  """
  Run a Hartree-Fock calculation with molsturm. The list of valid input
  parameters can be retrieved by the means of the list "hartree_fock_keys".

  Parameters:
    forward_parameters:   Should the list of kwargs be forwarded to the
                          output dictionary for archiving purposes?
                          If True(default) the returned dictionary
                          will contain a key "input_parameters" which
                          is a copy of the parameters on the commandline.

  A couple of selected kwargs:
    - basis_type        The type of the basis used for the calculation
    - coords            List of iterables of size 3: Coordinats of the atoms
    - atoms             List of the atom symbols (in the same order as coords)
  """

  # The list of valid keys is the list of keys
  # with the special ones (starting with __) removed.
  params_keys = [ k for k in dir(iface.Parameters) if k[0] != "_" ]
  res_keys = [ k for k in dir(iface.HfResults) if k[0] != "_" ]

  # These keys exist in the Parameters struct, but should not be exposed further up.
  exclude_keys = [ "restricted_set_by_user" ]

  # Build params and run:
  params = iface.Parameters()
  for key in kwargs:
    if not key in params_keys or key in exclude_keys:
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

