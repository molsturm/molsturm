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
import molsturm.posthf
import numpy as np


def ccd_residual(t2, fock, eri):
    """
    Compute the CCD residual tensor from the t2 amplitudes,
    the fock matrix and the repulsion integrals.
    """
    ret = np.einsum("abij->iajb", eri.block("vvoo")) \
        + np.einsum("ae,iejb->iajb", fock.block("vv"), t2) \
        - np.einsum("be,ieja->iajb", fock.block("vv"), t2) \
        - np.einsum("mi,majb->iajb", fock.block("oo"), t2) \
        + np.einsum("mj,maib->iajb", fock.block("oo"), t2)

    ret += \
        + 0.5 * np.einsum("mnij,manb->iajb", eri.block("oooo"), t2) \
        + 0.5 * np.einsum("abef,iejf->iajb", eri.block("vvvv"), t2) \
        + np.einsum("mbej,iame->iajb", eri.block("ovvo"), t2) \
        - np.einsum("mbei,jame->iajb", eri.block("ovvo"), t2) \
        - np.einsum("maej,ibme->iajb", eri.block("ovvo"), t2) \
        + np.einsum("maei,jbme->iajb", eri.block("ovvo"), t2)

    oovv = eri.block("oovv")
    ret += \
        - 0.5 * np.einsum("mnef,manf,iejb->iajb", oovv, t2, t2) \
        + 0.5 * np.einsum("mnef,mbnf,ieja->iajb", oovv, t2, t2) \
        - 0.5 * np.einsum("mnef,ienf,majb->iajb", oovv, t2, t2) \
        + 0.5 * np.einsum("mnef,jenf,maib->iajb", oovv, t2, t2) \
        + 0.25 * np.einsum("mnef,manb,iejf->iajb", oovv, t2, t2) \
        + 0.5 * np.einsum("mnef,iame,jbnf->iajb", oovv, t2, t2) \
        - 0.5 * np.einsum("mnef,jame,ibnf->iajb", oovv, t2, t2) \
        - 0.5 * np.einsum("mnef,ibme,janf->iajb", oovv, t2, t2) \
        + 0.5 * np.einsum("mnef,jbme,ianf->iajb", oovv, t2, t2)
    return ret


def ccd_approx_jacobian(t2, fock, eri):
    """
    Obtain an approximation for the CCD jacobian, i.e.
    for the derivative of the CCD residual wrt. t2 amplitudes.
    Note: This expression is only valid for canonical HF orbitals.
    """
    orben_o = fock.block("oo").diagonal()
    orben_v = fock.block("vv").diagonal()
    e_iajb = \
        - orben_o[:, np.newaxis, np.newaxis, np.newaxis] \
        + orben_v[np.newaxis, :, np.newaxis, np.newaxis] \
        - orben_o[np.newaxis, np.newaxis, :, np.newaxis] \
        + orben_v[np.newaxis, np.newaxis, np.newaxis, :]
    assert e_iajb.shape == t2.shape
    return e_iajb


def ccd_energy(t2, fock, eri):
    return 0.25 * np.einsum("ijab,iajb", eri.block("oovv"), t2)


def ccd(state):
    # Tolerance to which CCD is solved:
    tolerance = 100 * state["final_error_norm"]

    # Form antisymmetric repulsion tensor in Physicist's index convention:
    eri_phys = state.eri.transpose((0, 2, 1, 3))
    eri_asym = eri_phys - eri_phys.transpose((1, 0, 2, 3))

    # Compute MP2 t2 ampltudes as initial guess
    mp2 = molsturm.posthf.mp2(state)
    t2 = mp2["t2_ovov"]

    residual_norm = np.inf
    n_iter = 0
    print("    iter  CCD corr. energy   residual")
    fmt = "    {:2d}    {:15.12f}    {:6.3g}"
    while residual_norm > tolerance:
        # Compute current CCD energy and residual
        corr = ccd_energy(t2, state.fock, eri_asym)
        residual = ccd_residual(t2, state.fock, eri_asym)
        residual_norm = np.max(np.max(residual))
        print(fmt.format(n_iter, corr, residual_norm))

        # Quasi-Newton update step
        n_iter += 1
        t2 -= residual / ccd_approx_jacobian(t2, state.fock, eri_asym)

        if n_iter > 100:
            raise RuntimeError("CCD not converged after 100 iterations.")
    return corr, t2


if __name__ == "__main__":
    sys = molsturm.MolecularSystem(
        atoms=["O", "O"],
        coords=[(0, 0, 0), (0, 0, 2.8535)],
        multiplicity=3
    )
    state = molsturm.hartree_fock(sys, basis_type="gaussian",
                                  basis_set_name="6-31g", conv_tol=5e-7)

    corr, t2 = ccd(state)
    maxamp = np.max(np.abs(t2))
    print("CCD correlation energy: ", corr)
    print("Largest t2 amplitude:   ", maxamp)
    print()

    hf = state.energy_ground_state
    print("HF energy:              ", hf)
    print("CCD total energy:       ", corr + hf)
