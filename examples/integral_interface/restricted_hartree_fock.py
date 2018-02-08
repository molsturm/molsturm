#!/usr/bin/env python3

import molsturm
from molsturm import integrals
from scipy import linalg
import numpy as np

params = molsturm.ScfParameters()

# Supply molecular structure
atoms = ["Be"]  # as a list of names, symbols, atom numbers ...
positions = [[0, 0, 0]]  # as a list of triples x,y,z
system = molsturm.System(atoms, positions)
system.charge = 0  # Optional, 0 by default
system.multiplicity = 1  # Optional, 1 by default for even-electron systems,
#                                    2 for odd-electron systems
params.system = system

# Supply basis set information
basis = molsturm.construct_basis("gaussian", system, basis_set_name="pc-3")
# basis = molsturm.construct_basis("sturmian/atomic", system, k_exp=1.988,
#                                  n_max=5, l_max=1, m_max=1)
params.basis = basis

# Compute integral data
print("Computing integrals ... ", end="", flush=True)
s_bb = integrals.overlap_bb(params)
t_bb = integrals.kinetic_bb(params)
v_bb = integrals.nuclear_attraction_bb(params)
eri_bbbb = integrals.electron_repulsion_bbbb(params)
print("done")

#
# Restricted Hartree-Fock routine starts
#

# The number of occupied orbitals to compute
n_orb = system.n_alpha

# Form a hcore guess:
_, c_bf = linalg.eigh(t_bb + v_bb, s_bb, eigvals=(0, n_orb - 1))

oldene = 0
print("iter     etot          echange")
for i in range(1, 101):
    # Form new Fock matrix
    fock = t_bb + v_bb \
        + 2 * np.einsum("ai,bi,abcd->cd", c_bf, c_bf, eri_bbbb) \
        - 1 * np.einsum("ai,bi,cbad->cd", c_bf, c_bf, eri_bbbb)

    # Diagonalise it:
    enes, cf_bf = linalg.eigh(fock, s_bb)

    # Get occupied coefficients according to Aufbau principle
    c_bf = cf_bf[:, :n_orb]

    # Compute new HF energy
    ene_one_elec = np.einsum("ai,bi,ab", c_bf, c_bf, t_bb + v_bb)
    ene_two_elec = \
        2 * np.einsum("ai,bi,cj,dj,abcd", c_bf, c_bf, c_bf, c_bf, eri_bbbb) \
        - np.einsum("ai,bi,cj,dj,cbad", c_bf, c_bf, c_bf, c_bf, eri_bbbb)
    ene = 2 * (ene_one_elec + 0.5 * ene_two_elec)

    # Display current iteration
    ene_change = ene - oldene
    print(i, "  ", ene, "  ", ene_change)

    # Check for convergence
    if abs(ene_change) < 1e-6:
        break

    # Copy old energy
    oldene = ene
# end for
