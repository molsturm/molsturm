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


def distribute_electrons(n_electrons, multiplicity):
    """
    Return a tuple of (n_alpha, n_beta) resulting in the multiplicity
    and the electron count provided.
    If this cannot be achieved, raises a ValueError.
    """
    if multiplicity <= 0:
        raise ValueError("The multiplicity needs to be a positive number.")

    if n_electrons < 0:
        raise ValueError("The electron count needs to be zero or larger.")

    spin_twice = multiplicity - 1
    if spin_twice % 2 != n_electrons % 2:
        raise ValueError("Only a system with an even number of electrons can "
                         "have an odd multiplicity and vice versa. This system "
                         "has " + str(n_electrons) + " electrons, but a "
                         "multiplicity of " + str(multiplicity) + " is desired.")

    if spin_twice > n_electrons:
        raise ValueError("A system with " + str(n_electrons) + " electrons cannot "
                         "have a multiplicity larger than " + str(n_electrons + 1) +
                         ". You requested a multiplicity of " + str(multiplicity) +
                         ", however.")

    n_alpha = (n_electrons - spin_twice) // 2 + spin_twice
    n_beta = n_electrons - n_alpha

    assert n_alpha >= 0 and n_beta >= 0
    assert n_alpha - n_beta == spin_twice
    return n_alpha, n_beta


class System(Structure):
    """
    Class representing a molecular system, i.e. a structure and a number of electrons.

    The constructor of the class yields a neutral system with a multiplicity of
    1 (even-electron systems) or 2 (odd-electron systems).
    If this is not desired use the charge and multiplicity properties
    or the adjust_electrons method in order to yield a charged system or a triplet.

    For example. An O^+ quartet can made using this sequence of commands:
        >>> oxygen = System(["O"], [[0,0,0]])
        >>> oxygen.charge = 1
        >>> oxygen.multiplicity = 4

    Note that the order of the change and the multiplicity assignment is important,
    since changing the charge may change the multiplicity if for consistency
    this is needed.
    Similarly one could do the following

        >>> oxygen = System(["O"], [[0,0,0]])
        >>> oxygen.adjust_electrons(multiplicity=4, charge=1)

    Full flexibility exists, however, to make such adjustment by onself
    by providing the precise electron distribution upon construction.
        >>> oxygen = System(["O"], electrons=(5, 2))

    or in two steps by explicit assignment.
        >>> oxygen = System(["O"])
        >>> oxygen.n_alpha = 5
        >>> oxygen.n_beta = 2

    Related to this class is the method molsturm.read_system which reads
    a file, like an xyz file, and constructs an appropriate system object from it.
    """

    def __init__(self, atoms=[], coords=[], electrons=None):
        """
        Initialise a molecular system. By default it is empty (i.e. contains no atoms
        or electrons). In order to construct a non-empty system, specify one or more
        of the optional parameters.

        atoms:       List of all atoms of the structure. Either symbol or atomic numbers
                     or names of the elements are accepted.
        coords:      List of the coordinates of the involved atoms.
                     May be missing if none or only one atom is specified.
        electrons    The number of electrons in the system.
                     Can be given as a tuple (n_alpha, n_beta) or as a single
                     value. Then we take n_alpha = electrons // 2 and
                     n_beta = electrons - n_alpha
                     If absent, the value is set such that a neutral system
                     with multiplicity 1 (even electron systems)
                     or 2 (odd electron systems) results.

        In order to adjust the charge or multiplicity after construction,
        see the method adjust_electrons.
        """
        super().__init__(atoms, coords)
        if isinstance(electrons, (tuple, list)):
            self.n_alpha, self.n_beta = electrons
            if self.n_alpha < self.n_beta:
                raise ValueError("molsturm assumes that n_alpha >= n_beta at many "
                                 "places. You instructed to construct a system with "
                                 "n_beta (== " + str(self.n_beta) + ") > n_alpha ( == " +
                                 str(self.n_alpha) + "), which is not supported.")
            return

        if electrons is None:
            # Set electron count such that charge cancels to zero
            nuclear_charge = self.atom_numbers.sum()
            electrons = nuclear_charge
        elif not isinstance(electrons, int):
            raise TypeError("electrons needs to be a tuple, list or int")

        # Distribute to alpha and beta with alpha getting more electrons:
        self.n_beta = electrons // 2
        self.n_alpha = electrons - self.n_beta

        assert (electrons % 2 == 0 and self.multiplicity == 1) or \
            (electrons % 2 == 1 and self.multiplicity == 2)

    def adjust_electrons(self, charge=None, multiplicity=None,
                         allow_multiplity_change=True):
        """
        Distribute a certain number of electrons into alpha and beta electrons,
        such that a particular charge and multiplicity is achieved.

        charge         The total charge to achieve in the molecular system after the
                       distribution of electrons. If this parameter is absent,
                       the present value (by default zero) is not altered.
        multiplicity   The multiplicity to achieve after the distribution.
                       If this value is not set, the function will try to retain
                       the original multiplicity. If this cannot be done,
                       then the multiplicity will be 1 for even electron systems
                       and 2 for odd electron systems after the call.
        allow_multiplity_change   If this is set to false, than changing the charge
                       in a way that the current multiplicity value becomes
                       invalid is not allowed and yields a ValueError.
        """
        if charge is None:
            n_elec_count = self.n_electrons
        else:
            nuclear_charge = self.atom_numbers.sum()
            n_elec_count = nuclear_charge - charge
            if n_elec_count < 0:
                raise ValueError("Charge cannot be more positive than the "
                                 "total nuclear charge")

        if multiplicity is None:
            # Try first with original multiplicity.
            try:
                self.n_alpha, self.n_beta = distribute_electrons(n_elec_count,
                                                                 self.multiplicity)
                return
            except ValueError:
                if not allow_multiplity_change:
                    raise ValueError("Changing the charge in the way requested would "
                                     "require a change in multiplicity as well.")

            # Set to 1 for even-electron count and 2 for odd:
            multiplicity = 1 if n_elec_count % 2 == 0 else 2

        self.n_alpha, self.n_beta = distribute_electrons(n_elec_count, multiplicity)

    @property
    def multiplicity(self):
        """Return the multiplicity of the system"""
        return self.n_alpha - self.n_beta + 1

    @multiplicity.setter
    def multiplicity(self, value):
        self.adjust_electrons(multiplicity=value)

    @property
    def is_closed_shell(self):
        """Is the system closed shell"""
        return self.n_alpha == self.n_beta

    @property
    def n_electrons(self):
        return self.n_alpha + self.n_beta

    @n_electrons.setter
    def n_electrons(self, value):
        """Change the total number of electrons of the system.

        Setting this attribute as a single integer
        will change the multiplicity
        as well in case the original value cannot be retained.
        See adjust_electrons for details.
        """
        if isinstance(value, (tuple, list)):
            self.n_alpha, self.n_beta = value
        else:
            newcharge = self.atom_numbers.sum() - value
            self.adjust_electrons(charge=newcharge)

    @property
    def total_charge(self):
        """Return the total resulting charge of the system"""
        return self.atom_numbers.sum() - self.n_electrons

    @property
    def charge(self):
        """Return the total resulting charge of the system"""
        return self.total_charge

    @charge.setter
    def charge(self, value):
        """Change the total charge of the system.

        Setting this attribute will change the multiplicity
        as well in case the original value cannot be retained.

        See adjust_electrons for details.
        """
        self.adjust_electrons(charge=value)
