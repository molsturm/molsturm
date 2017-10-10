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

from ._scf import self_consistent_field
import numpy as np


def extrapolate_from_previous(old_state, scf_params):
    """
    Extrapolate the old SCF result state onto the new parameters to build
    a guess for a new scf procedure
    """
    scf_params = scf_params.copy()
    scf_params.normalise()
    old_params = old_state.input_parameters

    def check_agreement(key, message):
        if old_params[key] != scf_params[key]:
            raise ValueError(message + "(guess: " + str(old_params[key]) +
                             ", scfparams: " + str(scf_params[key]))

    check_agreement("scf/kind", "Cannot use restricted guess for unrestricted "
                    "calculation or vice versa.")
    check_agreement("discretisation/basis_type",
                    "discretisation/basis_type do not agree")

    old_orben_f = old_state["orben_f"]
    old_orbcoeff_bf = old_state["orbcoeff_bf"]
    try:
        # Project the old SCF result onto the new basis:
        P = old_params.basis.obtain_projection_to(scf_params.basis)
        proj_orbcoeff_bf = np.matmul(P, old_orbcoeff_bf)
    except (ValueError, NotImplementedError) as e:
        raise ValueError(str(e))

    # Check that sizes agree and build the actual guess:
    old_sizes = old_params.scf_sizes
    sizes = scf_params.scf_sizes
    assert old_sizes.n_spin == sizes.n_spin
    n_spin = sizes.n_spin

    assert old_orben_f.shape == (2 * old_sizes.n_fock, )
    assert proj_orbcoeff_bf.shape == (sizes.n_bas, 2 * old_sizes.n_fock)

    # Note: This truncation scheme works for both restricted and unrestricted
    #       cases due to the fact that the order is spin, bas, fock
    #       In other words the fock index runs fastest.

    # First reshape the old objects
    old_orben_f = old_orben_f[:old_sizes.n_fock * n_spin]
    proj_orbcoeff_bf = proj_orbcoeff_bf[:, :old_sizes.n_fock * n_spin]

    old_orben_f = old_orben_f.reshape(n_spin, old_sizes.n_fock)
    proj_orbcoeff_bf = proj_orbcoeff_bf.reshape(sizes.n_bas, n_spin, old_sizes.n_fock)
    proj_orbcoeff_bf = proj_orbcoeff_bf.transpose(1, 0, 2)

    # Then extend / truncate to the new expected shape
    orben_f = np.zeros((n_spin, sizes.n_fock))
    orbcoeff_bf = np.zeros((n_spin, sizes.n_bas, sizes.n_fock))

    min_fock = min(sizes.n_fock, old_sizes.n_fock)
    orben_f[:, :min_fock] = old_orben_f[:, :min_fock]
    orbcoeff_bf[:, :, :min_fock] = proj_orbcoeff_bf[:, :, :min_fock]

    return orben_f, orbcoeff_bf


def best_of_n(scfparams, n_repeats=4):
    """
    Find the best guess out of n_repeats random SCF calculations.

    Is intended to produce a good guess in case hcore and other
    methods fails
    """

    scfparams = scfparams.copy()
    scfparams["guess/method"] = "random"

    beststate = None
    for i in range(n_repeats):
        newstate = self_consistent_field(scfparams)

        if beststate is None or \
           newstate["energy_ground_state"] < beststate["energy_ground_state"]:
            beststate = newstate
    return extrapolate_from_previous(beststate, scfparams)
