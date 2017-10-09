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

        # Update geometry and run calculation:
        system = geometryfctn(args)
        res = molsturm.hartree_fock(system, conv_tol=conv_tol / 10,
                                    guess=recent_hf_res, **params)
        print("arguments: ", args, "niter:", res["n_iter"])

        # Update the stored recent_hf_res but only if it differs enough
        # from the one we currently have stored
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
    def construct_system(r_theta):
        r, theta = r_theta
        return molsturm.MolecularSystem(
            atoms=["O", "H", "H"],
            coords=[
                (0, 0, 0),
                (r, 0, 0),
                (r * np.cos(theta), r * np.sin(theta), 0)
            ]
        )
    return optimize_geometry(construct_system, (rHO_guess, angHO_guess), conv_tol,
                             level, **params)


if __name__ == "__main__":
    optimize_h2o(1.6, 105, conv_tol=1e-7, level="hf", basis_type="gaussian",
                 basis_set_name="def2-sv(p)")
