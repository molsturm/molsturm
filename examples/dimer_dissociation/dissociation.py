#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
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

import molsturm
import numpy as np


def compute_curve(atom, basis_set="sto-3g", conv_tol=1e-6, zrange=(0.5, 8.0), n_points=25,
                  restricted=False, verbose=False):
    z = np.linspace(zrange[0], zrange[1], n_points)
    f = np.empty_like(z)
    previous_hf = None

    for i in range(len(z) - 1, -1, -1):
        params = {
            "basis_type": "gaussian/libint",
            "basis_set":  basis_set,
            "atoms":      [atom, atom],
            "coords":     [(0, 0, 0), (0, 0, z[i])],
            "conv_tol":   conv_tol,
            "restricted": restricted,
            "print_iterations": verbose,
        }

        if previous_hf:
            params["guess"] = previous_hf
        else:
            params["guess"] = "random"

        try:
            hf = molsturm.hartree_fock(**params)
            f[i] = hf["energy_ground_state"]
            previous_hf = hf
        except RuntimeError as e:
            print("Caught error for z=", z)
            print(str(e))
            print()
            f[i] = np.nan

    return z, f
