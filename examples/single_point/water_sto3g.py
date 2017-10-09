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
from water import water


def run(**extra):
    res = molsturm.hartree_fock(water, basis_type="gaussian",
                                basis_set_name="sto-3g", **extra)

    molsturm.print_convergence_summary(res)
    molsturm.print_energies(res)
    molsturm.print_mo_occupation(res)
    molsturm.print_quote(res)
    return res


if __name__ == "__main__":
    run(print_iterations=True)
