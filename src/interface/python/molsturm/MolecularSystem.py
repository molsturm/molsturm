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

class MolecularSystem():
  """
  Class representing a molecular system
  """
  _fields = [ "atoms", "atom_numbers", "coords", "multiplicity", "charge" ]

  def __init__(self, atoms=None, atom_numbers=None, coords=None, multiplicity=None, charge=None):
    self.atoms = atoms
    self.atom_numbers = atom_numbers
    self.coords = coords
    self.multiplicity = multiplicity
    self.charge = charge

    if atoms is None and atom_numbers is None:
      raise ValueError("One of the parameters atoms or atom_numbers is required")
    if not atoms is None and not atom_numbers is None:
      raise ValueError("Only one of the parameters atoms or atom_numbers may be provided")


  def as_hartree_fock_parameters(self):
    """
    Return a dict which represents exactly this system as parameters
    for the hartree_fock function.
    """
    return { key : getattr(self, key) for key in self._fields
             if getattr(self, key) is not None }
