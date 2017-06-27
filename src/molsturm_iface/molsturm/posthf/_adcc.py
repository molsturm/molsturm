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

from .._constants import HFRES_ARRAY_KEYS

try:
  import adcc
  adcc_found = True
except ImportError:
  adcc_found = False

def generate_adcc_adcman_input(hfres):
  """
  Take the results dictionary from a hf calculation and build
  the input dictionary for a adcc.adcman calculation out of it.
  """
  include_keys = [ "n_alpha", "n_beta", "n_orbs_alpha", "n_orbs_beta",
                   "n_bas", "restricted", "threshold",
                   #
                   "orben_f", "eri_ffff", "fock_ff", "orbcoeff_fb",
                   #
                   "energy_nuclear_repulsion", "energy_nuclear_attraction",
                   "energy_coulomb", "energy_exchange", "energy_kinetic",
                 ]
  remap_keys = { "energy_total" : "energy_scf" }

  params = { k:hfres[k] for k in include_keys if k in hfres }
  params.update({ remap_keys[k]:hfres[k] for k in remap_keys
                  if k in hfres })
  return params

def run_adcc_adcman(hfres,**params):
  """
  Take hf results and extra parameters as an input
  and run adcc off it.
  """
  if not adcc_found:
    raise RuntimeError("Cannot run adcc.adcman: adcc not found.")

  # Check everything we need is there
  for k in HFRES_ARRAY_KEYS:
    if not k in hfres:
      raise ValueError("hfres parameters do not contain the required key '" +
                       k + "'.")

  for k in [ "method" ]:
    if not k in params:
      raise ValueError("The parameter key " + k + " is mandatory.")


  resp = generate_adcc_adcman_input(hfres)
  resp.update(params)
  return adcc.adcman(**resp)
