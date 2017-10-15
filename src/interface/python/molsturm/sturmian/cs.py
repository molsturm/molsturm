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

from .. import scf_guess
from .._scf import self_consistent_field
import numpy as np
import scipy.optimize
from ..MolecularSystem import MolecularSystem


def empirical_kopt(system):
    """
    Estimate the optimal value for k for Hartree-Fock
    based on an empirical formula and some known values.

    The result will be a good estimate, but still very far off.
    To obtain a better result, feed the obtained value for
    k into find_kopt.

    The input parameter may be a molsturm.MolecularSystem object
    or a molsturm.ScfParameters object.
    """
    if not isinstance(system, MolecularSystem):
        # Assume that we got an ScfParameters object
        system = system.system

    if len(system.atom_numbers) != 1:
        raise ValueError("Can only work on atomic systems")
    Z = system.atom_numbers[0]

    # TODO This is just something where we fitted some data to get this.
    #      I guess one could make some rationalisation here using the scaling of
    #      Z_eff and E_hydrogenic vs Z and such.
    # Intervals and the equations a*x + b
    intervals = [
        {"start": 1, "end":   2, "estimate": lambda n: 1.0},
        {"start": 2, "end":   3, "estimate": lambda n: 1.972},
        {"start": 3, "end":  11, "estimate": lambda n: 0.425 * n + 0.285},
        {"start": 11, "end":  19, "estimate": lambda n: 0.285 * n + 0.715},
        # Just out of my head
        {"start": 19, "end": 10000, "estimate": lambda n: 0.2 * n + 1.1},
    ]

    for interval in intervals:
        if interval["start"] <= Z and Z < interval["end"]:
            return float(interval["estimate"](Z))


def estimate_kopt(state):
    """
    Estimate an optimal k Coulomb sturmian exponent from
    a previously performed calculation.

    Currently only works well for HF calculations.

    The result will be a good estimate, but still very far off.
    To obtain a better result, feed the obtained value for
    k into find_kopt.
    """
    avgene = np.average(state.fock.block("oo").diagonal())
    k_hf = float(np.sqrt(-2 * avgene))
    return k_hf


def find_kopt(scfparams, conv_tol=1e-3, n_repeat_random=3, print_iterations=False,
              guess_update_threshold=0.1, optimisation_method="Brent",
              use_scfparams_exponent=False):
    """
    Find the optimal coulomb-sturmian exponent for Hartree-Fock.

    scfparams    Problem and basis description for which to
                 find the optimal sturmian exponent. An scf guess which
                 is specified in here will be taken into account.
                 Should be an ScfParameters object.
    conv_tol     Convergence tolerance for the optimal k
    n_repeat_random   Number of times solves with random initial guesses are
                      repeated in order to make the outcome more robust.
    print_iterations   Print messages regarding the iterations
    guess_update_threshold    Threshold in energy change below which
                              the guess will be updated from one
                              calcuation to another in the optimisation
                              procedure
    use_scfparams_exponent
         Use the exponent for the sturmians which is specified
         in the discretisation information of scfparams.

    Returns the state of the successful calculation or
    throws a RuntimeError if no convergence was achieved.
    """
    scfparams = scfparams.copy()

    if print_iterations:
        def opt_print(*msg):
            print(*msg)
    else:
        def opt_print(*msg):
            pass

    # If no guess for scf specified use a random guess:
    scfparams.setdefault("guess/method", "random")
    if scfparams.get("scf/max_iter", 100) <= 100:
        scfparams["scf/max_iter"] = 100

    # Do not solve too exactly in the beginning
    scfparams["scf/conv_tol"] = max(1e-5, conv_tol)

    #
    # Find initial guess values for k:
    #
    # Perform HF calculation on an estimate for k and fill
    # a dictionary with some important data
    def run_kjob(job):
        if "energy_ground_state" in job:
            return

        scfparams["discretisation/k_exp"] = job["k"]

        results = [self_consistent_field(scfparams)]
        if scfparams["guess/method"] == "random":
            for i in range(n_repeat_random - 1):
                # if the guess is random, do each calculation 3 times
                # to make sure we really do get the true SCF minimum
                results.append(self_consistent_field(scfparams))
            opt_print("Random energies for k=" + str(job["k"]) + ": ",
                      *[r["energy_ground_state"] for r in results])

        # Sort the results by energy
        results.sort(key=lambda x: x["energy_ground_state"])
        job["energy_ground_state"] = results[0]["energy_ground_state"]
        job["state"] = results[0]
        return job

    # Run both the guess encoded inside the scfparams and the empirical k_emp
    k_guess = scfparams["discretisation/k_exp"]
    k_emp = empirical_kopt(scfparams)

    if use_scfparams_exponent and k_guess != k_emp:
        job_guess = run_kjob({"k": k_guess})
        job_other = run_kjob({"k": k_emp})

        fmt_guess = "Using guesses {0:.3g} (from the scfparams) and " \
                    "{1:.3g} (from empirical_kopt)"
        opt_print(fmt_guess.format(k_guess, k_emp))
    else:
        job_guess = run_kjob({"k": k_emp})

        # Find another exponent by estimating from the HF results we just obtained
        k_hf = estimate_kopt(job_guess["state"])
        job_other = run_kjob({"k": k_hf})

        fmt_guess = "Using guesses {0:.3g} (from empirical_kopt) and " \
                    "{1:.3g} (from HF of the first)"
        opt_print(fmt_guess.format(k_emp, k_hf))

    # Build the search bracket:
    bracket = (job_guess["k"], job_other["k"])
    bracket = (min(bracket), max(bracket))

    # TODO If we have an actual guess, i.e. not need to resort to random stuff
    #      it will probably speed up the optimisation a lot if we find a proper
    #      bracket

    # Take the sate with the lowest energy as the guess
    guess = scf_guess.extrapolate_from_previous(
        min(job_guess, job_other, key=lambda x: x["energy_ground_state"])["state"],
        scfparams
    )
    scfparams.set_guess_external(*guess)

    hfres = None    # the HF result to return in the end
    hfguess = None  # the HF result to compute the guess from
    n_iter = 0

    # Objective function to optimise
    opt_print("  I       k n_scf_iter  energy_ground_state")
    fmt = "{0:3d} {1:7.5f}        {2:3d}    {3:12.8f}"

    def objective(x):
        nonlocal n_iter
        nonlocal hfres
        nonlocal hfguess
        nonlocal scfparams
        n_iter += 1

        scfparams["discretisation/k_exp"] = float(x)
        hfres = self_consistent_field(scfparams)
        opt_print(fmt.format(n_iter, x, hfres["n_iter"], hfres["energy_ground_state"]))

        # Update the guess, but only if the delta in energy is not too large
        # or if the delta in energy is not too small or too large
        if hfguess is None:
            diff = np.inf
        else:
            diff = abs(hfres["energy_ground_state"] - hfguess["energy_ground_state"])

            # Update conv_tol to gradually become more and more close to
            # what we want to achieve
            scfparams["scf/conv_tol"] = max(conv_tol / 1000, abs(diff))

        if diff < guess_update_threshold and diff > 10 * conv_tol:
            # Use the state which is the best state, i.e. which has
            # the lowest ground state energy
            if hfres["energy_ground_state"] < hfguess["energy_ground_state"]:
                hfguess = hfres
            guess = scf_guess.extrapolate_from_previous(hfguess, scfparams)
            scfparams.set_guess_external(*guess)
        return hfres["energy_ground_state"]

    oret = scipy.optimize.minimize_scalar(objective, method=optimisation_method,
                                          tol=conv_tol, bracket=bracket)

    # In older versions it's impossible to find out about the success
    if hasattr(oret, "success") and not oret.success:
        raise RuntimeError("Optimisation not successful: " + oret.message)
    else:
        if abs(oret.x - hfres.input_parameters.basis.k_exp) > conv_tol:
            scfparams["discretisation/k_exp"] = float(oret.x)
            hfres = self_consistent_field(scfparams)

    return hfres
