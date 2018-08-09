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

import numpy as np


# TODO More documentation:.
class soarray(np.ndarray):
    """Spin-Orbital array"""

    # Subclassing numpy arrays is a little different than usual.
    # See https://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    # for a guideline.

    def __new__(cls, n_orbs, n_electrons, inner):
        # Create numpy array instance
        order = "C" if not np.isfortran(inner) else "F"
        obj = super(soarray, cls).__new__(cls, shape=inner.shape, dtype=inner.dtype,
                                          buffer=inner, strides=inner.strides,
                                          order=order)

        if isinstance(n_orbs, int):
            n_orbs_beta = n_orbs // 2
            obj.n_orbs = (n_orbs - n_orbs_beta, n_orbs_beta)
        elif isinstance(n_orbs, tuple):
            obj.n_orbs = n_orbs
        else:
            raise TypeError("n_orbs needs to be a tuple or an int.")

        if isinstance(n_electrons, int):
            n_beta = n_electrons // 2
            obj.n_electrons = (n_electrons - n_beta, n_beta)
        elif isinstance(n_electrons, tuple):
            obj.n_electrons = n_electrons
        else:
            raise TypeError("n_electrons needs to be a tuple or an int.")

        if obj.n_orbs[0] < obj.n_electrons[0] or obj.n_orbs[1] < obj.n_electrons[1]:
            raise ValueError("Need less electrons than orbitals.")

        if any(sh != sum(obj.n_orbs) for sh in inner.shape):
            raise ValueError("Shape of the inner tensor needs to be equal to the number "
                             "of orbitals.")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        if not hasattr(obj, "n_electrons"):
            raise ValueError("Sorry, cannot construct soarray from a plain numpy ndarray")
        else:
            self.n_electrons = getattr(obj, "n_electrons")
            self.n_orbs = getattr(obj, "n_orbs")

    def block(self, block):
        if len(block) != self.ndim:
            raise ValueError("Tensor dimension and number of letters in "
                             "the block string needs to agree")

        n_alpha = self.n_electrons[0]
        n_beta = self.n_electrons[1]
        n_orbs_alpha = self.n_orbs[0]
        n_orbs_beta = self.n_orbs[1]

        occa = slice(0, n_alpha)
        occb = slice(n_orbs_alpha, n_orbs_alpha + n_beta)
        virta = slice(n_alpha, n_orbs_alpha)
        virtb = slice(n_orbs_alpha + n_beta, n_orbs_alpha + n_orbs_beta)

        ranges = {
            "o": (occa, occb),
            "v": (virta, virtb),
        }

        if any(a not in ranges for a in block):
            raise ValueError("Block string may only contain letters " +
                             ",".join(ranges.keys()))

        ret = self
        for i, b in enumerate(block):
            # Build numpy index:
            idxa = i * (slice(None), ) + (ranges[b][0], Ellipsis)
            idxb = i * (slice(None), ) + (ranges[b][1], Ellipsis)
            ret = np.concatenate((ret[idxa], ret[idxb]), axis=i)
        return ret
