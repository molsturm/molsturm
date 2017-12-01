#!/usr/bin/env python3
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
## vi: tabstop=2 shiftwidth=2 softtabstop=2 expandtab

import molsturm
import molsturm.posthf

sys = molsturm.System(["beryllium"], [[0, 0, 0]])
bas = molsturm.construct_basis("sturmian/atomic", sys, k_exp=2.1, n_max=11, l_max=0,
                               backend="cs_reference_pc")
res = molsturm.hartree_fock(sys, bas, conv_tol=1e-10, print_iterations=True)

molsturm.print_convergence_summary(res)
molsturm.print_energies(res)
molsturm.print_mo_occupation(res)
print()

res_adc = molsturm.posthf.mp2(res)
print("MP2 energy", res_adc["energy_mp2"])
print("tot energy", res_adc["energy_ground_state"])

molsturm.print_quote(res)
