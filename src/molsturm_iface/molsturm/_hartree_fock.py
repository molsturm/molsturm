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
from ._constants import HFRES_ARRAY_KEYS, INPUT_PARAMETER_KEY
from collections import Iterable
from .scf_guess import extrapolate_from_previous

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


def __setup_user_provided_guess(kwargs, inputargs):
  """
  This is triggered if guess_external_orben_f or guess_external_orbcoeff_bf
  is provided by the user. This may only be done if "guess" is also set
  to external. Alternatively we set guess to external and check that both
  values are provided.
  """
  if not "guess" in kwargs or kwargs["guess"] != "external":
    raise ValueError("The parameters 'guess_external_orben_f' or "
                     "'guess_external_orbcoeff_bf' may only be provided iff 'guess'"
                     " is 'external'")

  if not "guess_external_orben_f" in kwargs \
     or not "guess_external_orbcoeff_bf" in kwargs:
    raise ValueError("If the guess parameter is set to 'external', the parameters "
                     "'guess_external_orben_f' (for the guess orbital energies) and "
                     "'guess_external_orbcoeff_bf' (and the guess orthonormal SCF "
                     "coefficients) are both required.")

  inputargs["guess_external_orben_f"] = kwargs["guess_external_orben_f"]
  inputargs["guess_external_orbcoeff_bf"] = kwargs["guess_external_orbcoeff_bf"]
  inputargs["guess"] = "external"


def __setup_guess(kwargs, inputargs):
  """
  This is triggered if guess is provided by the user.
  """
  guess_key = "guess"
  guess = kwargs[guess_key]

  if isinstance(guess, str):
    if guess == "external":
      __setup_user_provided_guess(kwargs, inputargs)
    else:
      inputargs[guess_key] = guess
    return
  elif not isinstance(guess, dict):
    raise TypeError("The value provided to guess must be a previous hartree_fock result "+
                    " or a string describing a valid guess method.")

  inputargs["guess_external_orben_f"], inputargs["guess_external_orbcoeff_bf"] = \
      extrapolate_from_previous(guess, kwargs)
  inputargs["guess"] = "external"


# Map from input kwargs to inputargs list
__kwargs_parse_map = {
  "guess":                       __setup_guess,
  "guess_external_orben_f":      __setup_user_provided_guess,
  "guess_external_orbcoeff_bf":  __setup_user_provided_guess,
}

# Special input parameters for which the above conversion functions need
# to be used before assignment
__params_transform_map = {
  "coords":                      __to_coords,
  "atom_numbers":                __to_atom_numbers,
  "atoms":                       __to_atoms,
  "nlm_basis":                   __to_nlm,
  "guess_external_orben_f":      __to_double_vector,
  "guess_external_orbcoeff_bf":  __to_double_vector,
}


# Extra keys witch which we deal with on the python level
HF_EXTRA_KEYS = [ ]


# These keys exist in the Parameters struct, but should not be exposed
# via the hartree_fock interface
HF_PARAMS_EXCLUDE_KEYS = [ "all" ] + [ k for k in dir(iface.Parameters)
                                       if k[0] == "_" or k.startswith("internal_") ]

"""The list of keys understood by the hartree_fock function"""
hartree_fock_keys = HF_EXTRA_KEYS + [ k for k in dir(iface.Parameters)
                                      if not k in HF_PARAMS_EXCLUDE_KEYS ]

def hartree_fock(**kwargs):
  """
  Run a Hartree-Fock calculation with molsturm. The list of valid input
  parameters can be retrieved by the means of the list "hartree_fock_keys".

  A couple of selected kwargs:
    - basis_type        The type of the basis used for the calculation
    - coords            List of iterables of size 3: Coordinats of the atoms
    - atoms             List of the atom symbols (in the same order as coords)
    - guess             The guess method to use
  """
  #
  # Input
  #
  # Keys which need to be parsed *after* all other ones have been.
  delayed_keys = [ "guess", "guess_external_orben_f", "guess_external_orbcoeff_bf" ]

  inputargs = dict()
  for key in kwargs:
    if not key in hartree_fock_keys:
      raise ValueError("Keyword " + key + " is unknown to hartree_fock")
    if key in delayed_keys:
      continue

    # Copy key and (possibly transformed) value:
    if key in __kwargs_parse_map:
      __kwargs_parse_map[key](kwargs,inputargs)
    else:
      inputargs[key] = kwargs[key]

  # Deal with the delayed keys:
  for key in delayed_keys:
    if key in kwargs:
      __kwargs_parse_map[key](kwargs,inputargs)

  # Setup parameters:
  params = iface.Parameters()
  for key in inputargs:
    if key in __params_transform_map:
      setattr(params, key, __params_transform_map[key](inputargs[key]))
    else:
      setattr(params, key, inputargs[key])

  if "restricted" in kwargs:
    # Make a note that the user specified the restricted keyword
    params.internal_restricted_set_by_user = True

  #
  res = iface.hartree_fock(params)
  #

  #
  # Output
  #
  # TODO Better return a dict-like class instead of a dict. That way
  #      we can use the class more easily and distinguish between results
  #      at different levels better.
  res_keys = [ k for k in dir(iface.HfResults) if k[0] != "_" ]
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

  # Forward input parameters to output
  out[INPUT_PARAMETER_KEY] = inputargs

  return out


def compute_derived_hartree_fock_energies(hfres):
  """Compute various derived hartree-fock energy terms."""
  # TODO It would be better to have this in a hfres class,
  #      which is returned by the hartree_fock function
  res=dict()

  # Prefix all energy keys use:
  prefix = "energy_"

  # Classify the different keys:
  zeroElectron = [ "nuclear_repulsion" ]    # No electrons involved
  twoElectron = [ "coulomb", "exchange" ]   # 2 electron terms

  # Keys with special treatment
  special = zeroElectron + twoElectron + [ "ground_state" ]
  oneElectron = sorted([ k[len(prefix):] for k in hfres
                         if k.startswith(prefix) and \
                           not k[len(prefix):] in special
                       ])

  # All energy terms:
  energies = zeroElectron + oneElectron + twoElectron

  # Store individual terms in returned dictionary
  res["terms"] = { ene : hfres[prefix+ene] for ene in energies }

  # Derived energies:
  res[prefix + "ground_state"] = hfres[prefix + "ground_state"]
  res[prefix + "1e"]           = sum([ hfres[prefix+ene] for ene in oneElectron ])
  res[prefix + "2e"]           = sum([ hfres[prefix+ene] for ene in twoElectron ])
  res[prefix + "electronic"]   = res[prefix + "1e"] + res[prefix + "2e"]
  res[prefix + "nuclear"]      = hfres[prefix + "nuclear_repulsion"]
  res[prefix + "potential"]    = sum([ hfres[prefix+ene] for ene in energies
                                     if not ene in [ "kinetic" ] ])
  res[prefix + "kinetic"]      = hfres[prefix + "kinetic"]
  res["virial_ratio"]          = - res[prefix + "potential"] / res[prefix + "kinetic"]

  return res
