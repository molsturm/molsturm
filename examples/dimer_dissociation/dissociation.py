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
import molsturm.posthf
import numpy as np


def compute_curve(atom, basis_set="sto-3g", conv_tol=1e-6, zrange=(0.5, 8.0), n_points=25,
                  restricted=False, verbose=False, method="hf"):
    if method not in ["hf", "mp2"]:
        raise ValueError("Only implemented for hf and mp2")

    z = np.linspace(zrange[0], zrange[1], n_points)
    f = np.empty_like(z)
    previous_hf = None

    idcs = np.argsort(z)[::-1]
    for i in idcs:
        sys = molsturm.MolecularSystem(atoms=[atom, atom],
                                       coords=[(0, 0, 0), (0, 0, z[i])])

        guess = previous_hf if previous_hf is not None else "random"
        try:
            hf = molsturm.hartree_fock(sys, basis_type="gaussian", conv_tol=conv_tol,
                                       basis_set_name=basis_set, guess=guess,
                                       restricted=restricted, print_iterations=verbose,
                                       max_iter=100)

            if previous_hf is None:
                previous_hf = hf
            if method == "mp2":
                hf = molsturm.posthf.mp2(hf)

            f[i] = hf["energy_ground_state"]
        except RuntimeError as e:
            print("Caught error for z=", z[i])
            print(str(e))
            print()
            f[i] = np.nan

    return z, f
