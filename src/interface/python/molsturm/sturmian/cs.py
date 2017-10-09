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

from . import cs_kopt
from .._scf import self_consistent_field
from ..ScfParameters import ScfParameters


def find_kopt(scfparams=ScfParameters(), conv_tol=1e-3, keep_basis_k=False):
    """
    Find the optimal coulomb-sturmian exponent for Hartree-Fock.

    scfparams    Problem and basis description for which to
                 find the optimal sturmian exponent.
    conv_tol     Convergence tolerance for the optimal k
    keep_basis_k    Do not ignore the k_exp specified in scfparams.
                    Normally a crude guess is obtained from an
                    empirical formula.

    Returns the state of the successful calculation or
    throws a RuntimeError if no convergence was achieved.
    """
    scfparams = scfparams.copy()

    # Typically we need more iterations for these methods
    if scfparams.get("scf/max_iter", 100) <= 100:
        scfparams["scf/max_iter"] = 100

    if not keep_basis_k:
        k, _ = cs_kopt.empirical(scfparams)

    for method in ["hf_orben", "optimisation"]:
        scfparams["discretisation/k_exp"] = k
        k, ret = getattr(cs_kopt, method)(scfparams, conv_tol=conv_tol)

    if abs(k - ret.input_parameters.basis.k_exp) > conv_tol:
        scfparams["discretisation/k_exp"] = k
        ret = self_consistent_field(scfparams)

    return ret
