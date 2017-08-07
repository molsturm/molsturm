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
import scipy.optimize
import molsturm
import molsturm.posthf

def optimize_geometry(geometryfctn, geometry_params_guess, conv_tol,
                      level="hf", **params):
  # Start from hcore guess
  recent_hf_res = "hcore"
  recent_hf_res_for_args = None

  # The function to minimise
  def objective_function(args):
    nonlocal recent_hf_res
    nonlocal recent_hf_res_for_args

    hfparams = params
    hfparams["conv_tol"] = conv_tol / 10
    hfparams["atoms"], hfparams["coords"] = geometryfctn(args)
    hfparams["export_repulsion_integrals"] = level != "hf"

    res = molsturm.hartree_fock(**hfparams, guess=recent_hf_res)
    print("arguments: ", args, "niter:", res["n_iter"])

    # Update the hf_res but only if it differs enough
    if recent_hf_res_for_args is None or \
       np.max(np.abs(np.array(recent_hf_res_for_args) - np.array(args))) > conv_tol:
      recent_hf_res_for_args = args
      recent_hf_res = res
      print("  -> guess update")

    if level == "hf":
      return res["energy_ground_state"]
    elif level == "mp2":
      res = molsturm.posthf.mp2(res)
      return res["energy_ground_state"]
    else:
      raise NotImplementedError("Level of theory not yet implemented")

  res = scipy.optimize.minimize(objective_function, geometry_params_guess, tol=conv_tol)
  print(res)


def optimize_h2o(rHO_guess, angHO_guess, conv_tol, level="hf", **params):
  def geometry(r_theta):
    r, theta = r_theta
    return [ "O", "H", "H" ], [ (0,0,0), (r,0,0), (r*np.cos(theta), r*np.sin(theta), 0) ]
  return optimize_geometry(geometry, (rHO_guess, angHO_guess), conv_tol, level, **params)

if __name__ == "__main__":
  params = { "basis_type": "gaussian/libint", "basis_set": "def2-sv(p)" }
  optimize_h2o(1.6, 105, conv_tol=1e-7, level="hf", **params)
