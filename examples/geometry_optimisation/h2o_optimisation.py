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

method = "hf"
basis = "sto-3g"
params = {
  "atoms":         [ "O", "H", "H" ],
  "basis_type":    "gaussian/libint",
  "conv_tol":      1e-7,
}
recent_hf_res = None


def objective_function(r_theta):
  global recent_hf_res

  r, theta = r_theta
  print("r == ", r, " theta == ", theta)
  coords = [
    (0,0,0),
    (r,0,0),
    (r*np.cos(theta), r*np.sin(theta), 0),
  ]

  if recent_hf_res:
    guess = recent_hf_res
  else:
    guess = "hcore"

  if method == "hf":
    res = molsturm.hartree_fock(**params, coords=coords, guess=guess, basis_set=basis)
    return res["energy_ground_state"]
  elif method == "mp2":
    res = molsturm.hartree_fock(**params, export_repulsion_integrals=True,
                                coords=coords, guess=guess, basis_set=basis)
    res = molsturm.posthf.mp2(res)
    return res["energy_ground_state"]
  else:
    raise NotImplementedError("Level of theory not yet implemented")
  recent_hf_res = res

def run():
  res = scipy.optimize.minimize(objective_function, (1.3, 105), tol=1e-6)
  print(res)

if __name__ == "__main__":
  run()
