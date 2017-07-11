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

from .._constants import HFRES_OPTIONAL
from .._constants import HFRES_INPUT_PARAMETER_KEY
import numpy as np

try:
  import adcc
  adcc_found = True
except ImportError:
  adcc_found = False


# Keys to include verbatim from the hf res dictionary to the adcc input
__adcc_include_keys = [
  "n_alpha", "n_beta", "n_orbs_alpha", "n_orbs_beta",
  "n_bas", "restricted",
  #
  "orben_f", "eri_ffff", "fock_ff",
  #
  "energy_nuclear_repulsion", "energy_nuclear_attraction",
  "energy_coulomb", "energy_exchange", "energy_kinetic",
]

# Keys which are to be transformed and then included in the adcc dictionary
# First item of the tuple is the target key and the second is the
# transformation function.
__adcc_remap_keys = {
  "energy_ground_state" : ("energy_scf",  lambda x:  x                     ),
  #                                note: The copy is needed here, since the data
  #                                      needs to be in memory in contiguous stride
  "orbcoeff_bf"  : ("orbcoeff_fb", lambda x:  x.transpose().copy()  ),
}

# Parameters which are to be transformed
__adcc_remap_params = {
  "verbosity":   ("print_level", lambda x: x),
}


def generate_adcc_adcman_input(hfres):
  """
  Take the results dictionary from a hf calculation and build
  the input dictionary for a adcc.adcman calculation out of it.
  """
  try:
    params = { k:hfres[k] for k in __adcc_include_keys }

    # Remap both the key string as well as the value of all key-value pairs
    # in the __adcc_remap_keys dictionary.
    params.update({ __adcc_remap_keys[k][0] : __adcc_remap_keys[k][1](hfres[k])
                    for k in __adcc_remap_keys })

    params["threshold"] = 5*hfres["final_error_norm"]
  except KeyError as e:
    raise ValueError("The hartree_fock result dictionary does not contain the required " +
                     "key '" + e.args[0] + ".")

  return params

def run_adcc_adcman(hfres,**params):
  """
  Take hf results and extra parameters as an input
  and run adcc off it.
  """
  if not adcc_found:
    raise RuntimeError("Cannot run adcc.adcman: adcc not found.")

  # Check everything we need is there
  for key in [ key for key in HFRES_OPTIONAL
               if key in __adcc_include_keys or key in __adcc_remap_keys ]:
    if not key in hfres:
      raise ValueError("The hartree_fock result dictionary does not contain the required key '" +
                       key + "'. Please switch the parameter '" + HFRES_OPTIONAL[key] +
                       "' to True when running the molsturm.hartree_fock caculation.")

  for k in [ "method" ]:
    if not k in params:
      raise ValueError("The parameter " + k + " is required.")


  resp = generate_adcc_adcman_input(hfres)
  for p in params:
    if p in __adcc_remap_params:
      remap = __adcc_remap_params[p]
      resp[remap[0]] = remap[1](params[p])
    else:
      resp[p] = params[p]
  return adcc.adcman(**resp)

