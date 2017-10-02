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
import gint


system = molsturm.MolecularSystem(atoms=[4], coords=[[0, 0, 0]])
basis = gint.sturmian.atomic.Basis(system, k_exp=2.3, n_max=4, l_max=1)


def run_matrix_free(system, basis, extra={}):
    params = molsturm.ScfParameters()
    params.system = system
    params.basis = basis
    params["guess/method"] = "random"
    params["guess/eigensolver/method"] = "arpack"
    params["scf/eigensolver/method"] = "arpack"
    params["scf/n_eigenpairs"] = 2 * (basis.size // 2)
    params.update(extra)

    res = molsturm.self_consistent_field(params)
    molsturm.print_convergence_summary(res)
    molsturm.print_energies(res)
    molsturm.print_mo_occupation(res)
    molsturm.print_quote(res)
    return res


def run(basis, system):
    run_matrix_free(basis, system, {
        "scf/eigensolver/method": "lapack",
        "guess/eigensolver/method": "lapack"
    })


if __name__ == "__main__":
    run_matrix_free(system, basis)
