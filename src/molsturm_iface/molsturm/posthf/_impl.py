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

try:
  import pyadc
  _pyadc_found = True
except ImportError:
  _pyadc_found = False

def __build_available_methods():
  if _pyadc_found:
    return [ "mp2" ]
  else:
    return [ ]
available_methods = __build_available_methods()

def generate_pyadc_input(hfres):
  """
  Take the results dictionary from a hf calculation and build
  the input dictionary for a pyadc calculation out of it.
  """
  include_keys = [ "n_alpha", "n_beta", "n_orbs_alpha", "n_orbs_beta", 
                   "restricted", "threshold", 
                   #
                   "repulsion_integrals_ffff", "fock_ff", 
                   "coeff_fb", "orbital_energies_f", 
                   #
                   "energy_nuclear_repulsion", "energy_nuclear_attraction",
                   "energy_coulomb", "energy_exchange", "energy_kinetic",
                 ]
  remap_keys = { "energy_total" : "energy_scf" }

  params = { k:hfres[k] for k in include_keys if k in hfres }
  params.update({ remap_keys[k]:hfres[k] for k in remap_keys
                  if k in hfres })
  return params

def run_pyadc(hfres,**params):
  """
  Take hf results and extra parameters as an input
  and run pyadc off it.
  """
  if not _pyadc_found:
    raise RuntimeError("Cannot run pyadc: pyadc not found.")

  resp = generate_pyadc_input(hfres)
  resp.update(params)
  return pyadc.adc(**resp)

def mp2(hfres,**params):
  """
  Take hf results and extra parameters as an input
  and run an mp2 calculations. Return the resulting
  dictionary of computed data.
  """
  return run_pyadc(hfres, **params)

