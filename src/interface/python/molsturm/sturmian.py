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

import numpy as np

""" Some utility functions to deal with a sturmian basis """

class CoulombSturmianBasis:
  def __init__(self, k_exp, n_max, l_max=None, m_max=None, order="nlm"):
    """Build the sturmian basis from n_max, l_max, m_max.
       order: Use the specified ordering in the basis tuples.
    """

    if order != "nlm":
      raise NotImplementedError("Orders other than nlm are not implemented.")

    self.pure = True    # Only defined by n_max, l_max, m_max
    self.order = "nlm"

    if l_max is None: l_max = n_max-1
    if m_max is None: m_max = l_max
    self.n_max = n_max
    self.l_max = l_max
    self.m_max = m_max
    self.k_exp = k_exp

    self.functions = [
      (n,l,m) for n in range(n_max+1)
              for l in range(0,min(n,l_max+1))
              for m in range(-min(m_max,l),min(m_max,l)+1)
    ]


  def __len__(self):
    return len(self.functions)


  def __getitem__(self, x):
    return self.functions[x]


  def __eq__(self, other):
    return self.functions == other.functions


  def __ne__(self, other):
    return self.functions != other.functions


  def obtain_projection_to(self, other_basis):
    """Build a projection matrix, which interpolates this basis
       to another CoulombSturmianBasis object.

       I.e. if P is the result of this function, than
       np.matmul(P, old_orbcoeff_bf) computes the orbital
       coefficients in the new basis (== other_basis) from the
       ones in the old basis (== self)
    """
    if not isinstance(other_basis, CoulombSturmianBasis):
      raise TypeError("The other basis needs to be of type "
                      "CoulombSturmianBasis as well.")

    to_nlm = other_basis.functions
    from_nlm = self.functions
    if len(to_nlm) < len(from_nlm):
      raise NotImplementedError("The direction from larger to smaller "
                                "sturmian basis is not yet implemented.")

    proj = np.zeros((len(to_nlm), len(from_nlm)))
    for i in range(len(from_nlm)):
      bas_tuple = from_nlm[i]
      for j in range(len(to_nlm)):
        if to_nlm[j] == bas_tuple:
          proj[j, i] = 1
          break

    return proj


  def as_hartree_fock_parameters(self):
    """
    Return a dict which represents exactly this basis as parameters
    for the hartree_fock function.
    """
    if self.pure:
      copy_keys = [ "n_max", "l_max", "m_max", "k_exp" ]
      return { key : getattr(self, key) for key in copy_keys }
    else:
      return { "nlm_basis": self.functions, "k_exp": k_exp }
