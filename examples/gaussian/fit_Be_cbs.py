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
import molsturm.gaussian
import numpy as np
import matplotlib.pyplot as plt
import warnings


params = molsturm.ScfParameters()
params.system = molsturm.MolecularSystem("Be")
params["scf/print_iterations"] = True
params.basis = molsturm.construct_basis("gaussian", params.system,
                                        basis_set_name="cc-pvdz")

print("Running fit, please wait")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    res = molsturm.gaussian.extrapolate_cbs_limit(params)

# Print the result
cov = res.cov
print("Estimated parameters:  value (stddev)")
print("     E(CBS) == {0:16.14G} ({1: >9.4G})".format(res.cbs, np.sqrt(cov[0][0])))
print("          A == {0:16.14G} ({1: >9.4G})".format(res.A,   np.sqrt(cov[0][1])))
print("          B == {0:16.14G} ({1: >9.4G})".format(res.B,   np.sqrt(cov[1][2])))

# Plot the result
plt.plot(res.orders, res.quantities, 'kx', label="cc-pVnZ")
x = np.linspace(min(res.orders) - 0.1, 8, 200)
plt.plot(x, res.extrapolation_function(x), 'r-', label="Fitted convergence")
plt.legend()
plt.savefig("cbs_fit.pdf")
