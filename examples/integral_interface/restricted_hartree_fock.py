#!/usr/bin/env python3

import molsturm
import molsturm.sturmian
from molsturm import integrals
from scipy import linalg
import numpy as np

# Run a closed-shell hartree fock using a self-coded
# routine in python.


def run_scf(system, basis):
    params = molsturm.ScfParameters()
    params.system = system
    params.basis = basis

    if system.n_alpha != system.n_beta:
        raise ValueError("Sorry only closed shell")
    n_orb = system.n_alpha

    # get integrals in memory
    print("Computing integrals ... ", end="", flush=True)
    t_bb = integrals.kinetic_bb(params)
    v_bb = integrals.nuclear_attraction_bb(params)
    s_bb = integrals.overlap_bb(params)
    eri_bbbb = integrals.electron_repulsion_bbbb(params)
    print("done")

    # Fock update routine
    def update_fock(c_bf):
        return t_bb + v_bb \
            + 2 * np.einsum("ai,bi,abcd->cd", c_bf, c_bf, eri_bbbb) \
            - 1 * np.einsum("ai,bi,cbad->cd", c_bf, c_bf, eri_bbbb)

    def update_coefficients(fock):
        enes, cf_bf = linalg.eigh(fock, s_bb)
        c_bf = cf_bf[:, :n_orb]
        return c_bf

    def hf_energy(c_bf):
        one_elec = np.einsum("ai,bi,ab", c_bf, c_bf, t_bb + v_bb)
        two_elec = \
            2 * np.einsum("ai,bi,cj,dj,abcd", c_bf, c_bf, c_bf, c_bf, eri_bbbb) \
            - np.einsum("ai,bi,cj,dj,cbad", c_bf, c_bf, c_bf, c_bf, eri_bbbb)
        return 2 * (one_elec + 0.5 * two_elec)

    def form_guess():
        # Get an hcore guess:
        _, c_bf = linalg.eigh(t_bb + v_bb, s_bb, eigvals=(0, n_orb - 1))
        assert c_bf.shape == (basis.size, n_orb)
        return c_bf

    c_bf = form_guess()
    oldene = 0
    print("iter     etot          echange")
    for i in range(1, 101):
        fock = update_fock(c_bf)
        c_bf = update_coefficients(fock)

        # Convergence display and check
        ene = hf_energy(c_bf)
        change = ene - oldene
        print(i, "  ", ene, "  ", change)
        oldene = ene
        if abs(change) < 1e-6:
            break


def main():
    system = molsturm.System(["Be"])
    kopt = molsturm.sturmian.cs.empirical_kopt(system)
    basis = molsturm.construct_basis("gaussian", system, basis_set_name="pc-3")
    run_scf(system, basis)


if __name__ == "__main__":
    main()
