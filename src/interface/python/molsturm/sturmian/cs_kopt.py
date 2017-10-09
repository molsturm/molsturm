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
import scipy.optimize

"""
This file contains routines for finding the optimal
k exponent of a coulomb-sturmian basis set
"""


def optimisation(scfparams, conv_tol=1e-3,
                 method="Brent", guess_previous_maxdelta=0.1):
    """
    Perform a full optimisation on the E vs k surface.

    scfparams  The problem as well as the discretisation options
               for which to find the optimal k
    conv_tol   The convergence tolerance for k
    method     The scalar optimisation method to use
    guess_previous_maxdelta    The maximal difference in energy
                               in order for a particular hf calculation
                               to be used as a guess for the next.

    Returns the predicted optimal k as well as the most recent
    calculation result which lead to it
    """
    scfparams = scfparams.copy()

    hfres = None  # the HF result to return in the end
    guess = None  # The guess to use
    k_guess = scfparams["discretisation/k_exp"]

    def objective(x):
        nonlocal hfres
        nonlocal guess
        nonlocal scfparams

        if guess is None:
            scfparams["guess/method"] = "random"

            # Do not solve too exactly in the begining
            scfparams["scf/conv_tol"] = max(1e-6, conv_tol)
        else:
            scfparams.set_guess_external(*guess)
            scfparams["scf/conv_tol"] = conv_tol / 100

        scfparams["discretisation/k_exp"] = float(x + k_guess)
        hfnew = molsturm.self_consistent_field(scfparams)

        # Update the guess, but only if the delta in energy is not too large
        # or if the delta in energy is not too small
        if hfres is not None:
            delta = abs(hfnew["energy_ground_state"] - hfres["energy_ground_state"])
        else:
            delta = 10 * guess_previous_maxdelta
        if delta < guess_previous_maxdelta and delta > 10 * conv_tol:
            guess = molsturm.scf_guess.extrapolate_from_previous(hfnew, scfparams)
        hfres = hfnew
        return hfres["energy_ground_state"]

    oret = scipy.optimize.minimize_scalar(objective, method=method, tol=conv_tol)
    if not oret.success:
        raise RuntimeError("Optimisation not successful: " + oret.message)
    else:
        return float(oret.x + k_guess), hfres


def hf_orben(scfparams, conv_tol=1e-3):
    """
    Run a HF calculation in order to determine the next k according to
    the formula E = -0.5*k^2, where E is the average HF orbital energy.
    Even though this method is extremely simple, it works remarkably well
    for getting into the correct range for k and should hence be used to
    improve a first guess. Usually one or two sweeps are enough as the
    effect quickly diminishes.

    Returns the predicted optimal k as well as the result of the
    calculation which lead to it.
    """

    scfparams = scfparams.copy()
    scfparams["guess/method"] = "random"
    scfparams["scf/conv_tol"] = conv_tol

    hfres = molsturm.self_consistent_field(scfparams)
    avgene = np.average(hfres.fock.block("oo").diagonal())
    k = float(np.sqrt(-2 * avgene))

    return k, hfres


def empirical(scfparams):
    """
    Obtain an optimal value for k_exp by looking up in a table
    or by considering an empirical formula for the atomic number.

    Returns the predicted optimal k as well as None.
    """
    system = scfparams.system
    assert len(system.atom_numbers) == 1
    Z = system.atom_numbers[0]

    # TODO This is just something out of my little hat
    #      I guess one could make some rationalisation here
    #      using the scaling of Z_eff and E_hydrogenic vs Z
    #      and such.
    k = float(np.sqrt(Z))

    return k, None
