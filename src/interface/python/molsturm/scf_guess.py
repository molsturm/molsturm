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
import warnings


def extrapolate_from_previous(old_state, scf_params):
    """
    Extrapolate the old SCF result state onto the new parameters to build
    a guess for a new scf procedure
    """
    warnings.warn("External guess of guess from previous is not yet "
                  "re-implemented properly.")

    # TODO This is just to make it work ... we need much more checking here.
    #      See the old version in scf_guess.old.py for ideas.

    scf_params = scf_params.copy()
    scf_params.normalise()

    scf_sizes = scf_params.scf_sizes
    n_spin = scf_sizes.n_spin
    n_fock = scf_sizes.n_fock
    n_bas = scf_sizes.n_bas

    if old_state["restricted"]:
        if not n_spin == 1:
            raise ValueError("Restricted guess for unrestricted computation")

        orben_f = old_state["orben_f"][:n_fock].reshape(n_spin, n_fock)
        orbcoeff_bf = old_state["orbcoeff_bf"][:, :n_fock]
        orbcoeff_bf = orbcoeff_bf.reshape(n_spin, n_bas, n_fock)
    else:
        if not n_spin == 2:
            raise ValueError("Unrestricted guess for restricted computation")

        orben_f = old_state["orben_f"].reshape(n_spin, n_fock)

        orbcoeff_bf = old_state["orbcoeff_bf"].reshape(n_bas, n_spin, n_fock)
        orbcoeff_bf = orbcoeff_bf.transpose(1, 0, 2)

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
