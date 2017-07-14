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
from .._constants import INPUT_PARAMETER_KEY
from .._basis import has_real_harmonics

try:
  import pyscf
  pyscf_found = True
except ImportError:
  pyscf_found = False

# Remap for parameters
__pyscf_remap_fci_keys = {
  "n_roots":   ("nroots", lambda x:x),
  "verbosity": ("verbose", lambda x:x),
}

def run_fci_pyscf(hfres, **params):
  import pyscf.fci

  required_keys = [ "eri_ffff", "hcore_ff" ]

  rhf = hfres["restricted"]
  n_elec = hfres["n_alpha"] + hfres["n_beta"]

  for key in required_keys:
    if not key in hfres:
      raise ValueError("The hartree_fock result dictionary does not contain the required key '" +
                       key + "'. Please switch the parameter '" + HFRES_OPTIONAL[key] +
                       "' to True when running the molsturm.hartree_fock caculation.")


  # If the underlying basis uses complex harmonics than we cannot
  # assume that (ij|kl) = (ji|kl) (Shell-pair, Mullikan notation)
  real_harmonics = has_real_harmonics(**hfres[INPUT_PARAMETER_KEY])

  if rhf and real_harmonics:
    fci = pyscf.fci.direct_spin1.FCI()
    spin_square_multipl = pyscf.fci.spin_op.spin_square
    #
    # Only use the alpha-alpha block
    n_orbs = hfres["n_orbs_alpha"]
    h1e = hfres["hcore_ff"][:n_orbs,:n_orbs]
    eri = hfres["eri_ffff"][:n_orbs,:n_orbs,:n_orbs,:n_orbs]
  elif not rhf and real_harmonics:
    fci = pyscf.fci.direct_uhf.FCI()
    spin_square_multipl = pyscf.fci.spin_op.spin_square
    #
    # Use alpha and beta blocks
    n_orbs = hfres["n_orbs_alpha"] + hfres["n_orbs_beta"]
    h1e = hfres["hcore_ff"]
    eri = hfres["eri_ffff"]
    raise NotImplementedError("UHF not yet implemented.")
  elif rhf and not real_harmonics:
    # The only solver pyscf offers for this case is the
    # generic nosym solver where (ij|kl) = (ji|kl) is not
    # assumed
    fci = pyscf.fci.direct_nosym.FCI()
    spin_square_multipl = pyscf.fci.spin_op.spin_square
    #
    # Only use the alpha-alpha block
    n_orbs = hfres["n_orbs_alpha"]
    h1e = hfres["hcore_ff"][:n_orbs,:n_orbs]
    eri = hfres["eri_ffff"][:n_orbs,:n_orbs,:n_orbs,:n_orbs]
  elif not rhf and not real_harmonics:
    raise ValueError("Currently doing FCI for an UHF ground state "
                     "and without real spherical harmonics "
                     "(Sturmians) is not supported.")
  else:
    assert False, "Internal error"

  # Set the parameters:
  for p in params:
    if p in __pyscf_remap_fci_keys:
      remap = __pyscf_remap_fci_keys[p]
      setattr(fci, remap[0], remap[1](params[p]))
    else:
      setattr(fci, p, params[p])

  if not "conv_tol" in params:
    setattr(fci, "conv_tol", 5*hfres["final_error_norm"])

  result = fci.kernel(h1e, eri, n_orbs, (hfres["n_alpha"], hfres["n_beta"]),
                      ecore=hfres["energy_nuclear_repulsion"])

  # result is a tuple (energy, fcivector) iff nroots == 1,
  # else it is a list of these objects
  #   => Normalise here
  if params.get("n_roots", 1) == 1:
    result = [ result[0] ], [ result[1] ]

  # TODO I do not quite understand the ci vector object
  #      and what it means at the moment.

  res={ "states": [], INPUT_PARAMETER_KEY: params }
  for i in range(len(result[0])):
    civector = result[1][i]
    state = {
      "energy":   result[0][i],
      "civector": civector,
    }
    state["spin_squared"], state["multiplicity"] = \
        spin_square_multipl(civector, n_orbs, n_elec)
    res["states"].append(state)

  return res

