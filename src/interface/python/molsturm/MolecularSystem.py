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

from gint import Structure


class MolecularSystem(Structure):
    """
    Class representing a molecular system, i.e. a structure and a number of electrons.
    """
    def __init__(self, atoms=None, coords=None, electrons=None, multiplicity=None,
                 charge=None, structure=None):
        """
        Initialise a molecular system. By default it is empty (i.e. contains no atoms
        or electrons). In order to construct a non-empty system, specify one or more
        of the optional parameters. Certain combinations are not valid and will
        raise a ValueError or TypeError.

        electrons       The number of electrons in the system.
                        Can be given as a tuple (n_alpha, n_beta) or as a single
                        value. Then we take n_alpha = electrons // 2 and
                        n_beta = electrons - n_alpha
        multiplicity    The multiplicity (2S+1) of the ground state to compute.
                        By default 1 is chosen for even electron systems and 2
                        for odd electron systems.
        charge          Computed from the electrons and the nuclei in the system
                        or set to 0 (neutral atom).
        """

        if structure is None:
            super().__init__(atoms, coords)
        elif not isinstance(structure, Structure):
            raise TypeError("structure needs to be a gint.Structure object.")
        else:
            super().__init__(structure.atom_numbers, structure.coords)

        if electrons is None:
            if charge is None:
                charge = 0
            n_elec_count = super().total_charge - charge

            if multiplicity is None:
                if n_elec_count % 2 == 0:
                    multiplicity = 1
                else:
                    multiplicity = 2
            elif multiplicity <= 0:
                raise ValueError("The multiplicity needs to be a positive number.")
            spin_twice = multiplicity - 1

            if spin_twice % 2 != n_elec_count % 2:
                raise ValueError("Only a system with an even number of electrons can "
                                 "have an odd multiplicity and vice versa. This system "
                                 "has " + str(n_elec_count) + " electrons, but a "
                                 "multiplicity of " + str(multiplicity) + ".")

            self.n_alpha = (n_elec_count - spin_twice) // 2 + spin_twice
            self.n_beta = n_elec_count - self.n_alpha

            assert self.multiplicity == multiplicity
            assert self.n_alpha - self.n_beta == spin_twice
        elif isinstance(electrons, tuple):
            self.n_alpha, self.n_beta = electrons
        else:
            self.n_alpha = electrons // 2
            self.n_beta = electrons - self.n_alpha

        self.__charge = super().total_charge - self.n_alpha - self.n_beta

    @property
    def multiplicity(self):
        """Return the multiplicity of the system"""
        return self.n_alpha - self.n_beta + 1

    @property
    def is_closed_shell(self):
        """Is the system closed shell"""
        return self.n_alpha == self.n_beta

    @property
    def total_charge(self):
        """Return the total resulting charge of the system"""
        return self.__charge

    @property
    def charge(self):
        return self.__charge
    charge.__doc__ = total_charge.__doc__
