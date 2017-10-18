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
    scf_params.clear_guess()
    scf_params.normalise()
    old_params = old_state.input_parameters

    def check_agreement(key, message):
        if old_params[key] != scf_params[key]:
            raise ValueError(message + "(guess: " + str(old_params[key]) +
                             ", scfparams: " + str(scf_params[key]))

    if old_params["scf/kind"] != scf_params["scf/kind"]:
        raise ValueError("Cannot use restricted guess for unrestricted " +
                         "calculation or vice versa." +
                         "(guess scf/kind: " + str(old_params["scf/kind"]) +
                         ", scfparams: " + str(scf_params["scf/kind"]))

    if not isinstance(old_params.basis, type(scf_params.basis)):
        raise ValueError("The type of the basis used in guess and scfparams " +
                         "has to agree. Note that the precise backend may differ." +
                         "guess discretisation/basis_type: " +
                         str(old_params["discretisation/basis_type"]) +
                         ", scfparams: " +
                         str(scf_params["discretisation/basis_type"]))

    old_orben_f = old_state["orben_f"]
    old_orbcoeff_bf = old_state["orbcoeff_bf"]
    try:
        # Project the old SCF result onto the new basis:
        P = old_params.basis.obtain_projection_to(scf_params.basis)
        proj_orbcoeff_bf = np.dot(P, old_orbcoeff_bf)
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


def best_of_n(scfparams, n_repeats=4, n_max_tries=None):
    """
    Find the best guess out of n_repeats random SCF calculations.

    Is intended to produce a good guess in case hcore and other
    methods fails.

    n_repeats   Number of scf minima to compute before selecting the
                lowest-energy one
    n_max_tries   Number of tries to attempt. for getting an scf minimum.
                  After this number the procedure fails.

    If a run fails (e.g. because of max_iter is reached) then
    it is discarded. After n_max_tries without n_repeats successful
    SCFs, the most recent error is raised.
    """

    scfparams = scfparams.copy()
    scfparams.clear_guess()
    scfparams["guess/method"] = "random"

    # Set to at least 100 iterations
    scfparams.setdefault("scf/max_iter", 100)
    scfparams["scf/max_iter"] = max(scfparams["scf/max_iter"], 100)

    if n_max_tries is None:
        n_max_tries = 2 * n_repeats

    # Since the scf runs could fail (we are using a random guess after all)
    # We actually run up to n_max_tries times and skip erroneous runs

    beststate = None
    n_successful = 0
    for itry in range(1, n_max_tries + 1):
        try:
            newstate = self_consistent_field(scfparams)
            n_successful += 1
        except RuntimeError as e:
            if itry >= n_max_tries:
                raise e
            else:
                continue

        if beststate is None or \
           newstate["energy_ground_state"] < beststate["energy_ground_state"]:
            beststate = newstate

        if n_successful == n_repeats:
            break

    assert n_successful == n_repeats
    return extrapolate_from_previous(beststate, scfparams)
