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

def build_nlm_basis(n_max, l_max=None, m_max=None, order="nlm"):
  """Build the (n,l,m) tuples for a basis with a
     given n_max, l_max, n_max

     order: Use the specified ordering in the basis tuples.
  """
  if order != "nlm":
    raise NotImplementedError("Only order == 'nlm' is implemented.")

  if l_max is None: l_max = n_max-1
  if m_max is None: m_max = l_max
  return [ (n,l,m) for n in range(n_max+1)
                   for l in range(0,min(n,l_max+1))
                   for m in range(-min(m_max,l),min(m_max,l)+1)
         ]

def build_basis_projector(from_nlm,to_nlm):
  """
  Build a projector which projects from one nlm tuple list
  of Coulomb-Sturmian basis functions to another.

  I.e. if P is the result of this
  function, than np.matmul(P,old_orbcoeff_bf) computes the orbital
  coefficients in the new basis (== to_nlm) from the ones in the old
  basis (== from_nlm)
  """
  if len(to_nlm) < len(from_nlm):
    raise NotImplementedError("The direction from larger to smaller "
                              "sturmian basis is not yet implemented.")

  proj = np.zeros((len(to_nlm), len(from_nlm)))
  for i in range(len(from_nlm)):
    bas_tuple = from_nlm[i]
    for j in range(len(to_nlm)):
      if to_nlm[j] == bas_tuple:
        proj[j,i] = 1
        break

  return proj

