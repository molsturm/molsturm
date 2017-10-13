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

from ._constants import INPUT_PARAMETER_KEY
from .soarray import soarray
from .ScfParameters import ScfParameters

# TODO See also the ideas in StateNew.py

# TODO Quick and dirty wrapper to mimic the future syntax
class TmpState(dict):
    def __init__(self, other):
        for key in other:
            self.__setitem__(key, other[key])

    def eri_block(self, block="ffff"):
        """Get a block of the electron repulsion tensor
           in shell pair or chemists notation.
           No antisymmetrisation done yet.
        """
        if block != "ffff":
            raise NotImplementedError("Only ffff can be obtained at the moment.")
        return self.make_soarray(self.__getitem__("eri_ffff"))

    @property
    def eri(self):
        """
        Return the full eri tensor.
        """
        return self.eri_block()

    @property
    def fock(self):
        """
        Return the full fock matrix.
        """
        return self.make_soarray(self.__getitem__("fock_ff"))

    @property
    def energy_ground_state(self):
        return self.__getitem__("energy_ground_state")

    @property
    def n_orbs(self):
        return (self.__getitem__("n_orbs_alpha"), self.__getitem__("n_orbs_beta"))

    @property
    def n_electrons(self):
        return (self.__getitem__("n_alpha"), self.__getitem__("n_beta"))

    @property
    def input_parameters(self):
        return ScfParameters.from_dict(self[INPUT_PARAMETER_KEY])

    def make_soarray(self, inner):
        a = soarray(self.n_orbs, self.n_electrons, inner)
        return a


State = TmpState
