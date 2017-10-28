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

# This file contains the high-level fronted for SCF calculations
# which is around for convenience

from ._scf import self_consistent_field
from .scf_guess import extrapolate_from_previous
from .ScfParameters import ScfParameters
from .State import State
import gint
import numpy as np


def hartree_fock(system, basis=None, basis_type=None,
                 conv_tol=None, max_iter=None, n_eigenpairs=None,
                 restricted=None, guess=None, eigensolver=None,
                 print_iterations=None, **kwargs):
    """
    Run a Hartree-Fock calculation with molsturm.

    system        The molecular system to model
    basis         A valid basis object. If None the basis will be constructed
                  on the fly from teh basis_type and the kwargs.
    basis_type    String describing the type of basis function to use.
    conv_tol      SCF convergence tolerance
    max_iter      Maximum number of SCF iterations
    n_eigenpairs  Number of orbitals to compute
    eigensolver   SCF eigensolver to use
    print_iterations  Shall some diagnostics about the SCF iterations be printed
    restricted    Shall an restricted or an unrestricted SCF be run
    guess         Use another state or chose the guess method (as a string)

    Examples:
        hartree_fock(system="Be", basis_type="gaussian",
                     basis_set_name="sto-3g")

    For a more low-level entry point, which offers a larger range of parameters to
    influence the bahaviour of an SCF, see the function self_consistent_field.
    This function "hartree_fock" is essentially a wrapper around
    self_consistent_field, which performs roughly the steps:
        basis = molsturm.construct_basis(system, **kwargs)
        params = ScfParameters.from_args(system, basis, conv_tol, max_iter, ...)
        self_consistent_field(params)
    """
    if basis is None:
        if basis_type is None:
            raise ValueError("Either the basis or the basis_type needs to be given.")
        basis = gint.construct_basis(basis_type, system, **kwargs)

    # Construct guess from the HF arguments
    params = ScfParameters.from_args(system, basis, conv_tol, max_iter, n_eigenpairs,
                                     restricted, eigensolver, print_iterations)
    params.normalise()

    # Build and add guess parameters
    if isinstance(guess, str):
        if guess == "external":
            raise ValueError("External guess is only available via the "
                             "ScfParameters object.")
        else:
            params["guess/method"] = str(guess)
    elif isinstance(guess, State):
        # Normalise parameter tree such that the extrapolate_from_previous
        # has access to a proper set of parameters.
        params.normalise()

        # Extrapolate and add the external guess:
        orben_f, orbcoeff_bf = extrapolate_from_previous(guess, params)
        params.set_guess_external(orben_f, orbcoeff_bf)
    elif guess is None:
        pass
    else:
        raise TypeError("guess can only be a string or a State object.")

    return self_consistent_field(params)


def compute_derived_hartree_fock_energies(hfres):
    """Compute various derived hartree-fock energy terms."""
    # TODO It would be better to have this in a hfres class,
    #      which is returned by the hartree_fock function
    res = dict()

    # Prefix all energy keys use:
    prefix = "energy_"

    # Classify the different keys:
    zeroElectron = ["nuclear_repulsion"]    # No electrons involved
    twoElectron = ["coulomb", "exchange"]   # 2 electron terms

    # Keys with special treatment
    special = zeroElectron + twoElectron + ["ground_state"]
    oneElectron = sorted([k[len(prefix):]
                          for k in hfres
                          if k.startswith(prefix) and not k[len(prefix):] in special
                          ])

    # All energy terms:
    energies = zeroElectron + oneElectron + twoElectron

    # Store individual terms in returned dictionary
    res["terms"] = {ene: hfres[prefix + ene] for ene in energies}

    # Derived energies:
    res[prefix + "ground_state"] = hfres[prefix + "ground_state"]
    res[prefix + "1e"] = sum([hfres[prefix + ene] for ene in oneElectron])
    res[prefix + "2e"] = sum([hfres[prefix + ene] for ene in twoElectron])
    res[prefix + "electronic"] = res[prefix + "1e"] + res[prefix + "2e"]
    res[prefix + "nuclear"] = hfres[prefix + "nuclear_repulsion"]
    res[prefix + "potential"] = sum([hfres[prefix + ene] for ene in energies
                                     if ene not in ["kinetic"]])
    res[prefix + "kinetic"] = hfres[prefix + "kinetic"]
    res["virial_ratio"] = - res[prefix + "potential"] / res[prefix + "kinetic"]

    return res


def compute_coulomb_ff(hfres):
    """Compute the coulomb matrix in MO space"""
    noa = hfres["n_orbs_alpha"]
    na = hfres["n_alpha"]
    nb = hfres["n_beta"]
    jirep = hfres["eri_ffff"]

    return np.trace(jirep[:na, :na, :, :], axis1=0, axis2=1) + \
        np.trace(jirep[noa:noa + nb, noa:noa + nb, :, :], axis1=0, axis2=1)


def compute_exchange_ff(hfres):
    """Compute the exchange matrix in MO space"""
    noa = hfres["n_orbs_alpha"]
    na = hfres["n_alpha"]
    nb = hfres["n_beta"]
    jirep = hfres["eri_ffff"]

    return np.trace(jirep[:na, :, :, :na], axis1=0, axis2=3) + \
        np.trace(jirep[noa:noa + nb, :, :, noa:noa + nb], axis1=0, axis2=3)
