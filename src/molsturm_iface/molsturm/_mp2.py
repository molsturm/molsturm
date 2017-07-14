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

from ._utils import occ_range, virt_range
import numpy as np

def mp2(hfres, **params):
  """
  Compute the mp2 energy making no assumptions about the repulsion
  integrals (i.e. this works for both restricted as well as unrestricted)

  Similarly no assumptions about the structure of the repulsion integrals
  are made.

  Returns the mp2 energy as well as the computed t2 amplitude tensor.
  """
  # Setup an empty T2 amplitude tensor
  n_orbs = hfres["n_orbs_alpha"] + hfres["n_orbs_beta"]
  t2 = np.zeros((n_orbs, n_orbs, n_orbs, n_orbs))

  emp2 = 0  # Final mp2 energy contribution

  orben = hfres["orben_f"]
  eri = hfres["eri_ffff"]

  # TODO
  # This is the least efficient way possible:
  for i in occ_range(hfres):
    for j in occ_range(hfres):
      for a in virt_range(hfres):
        for b in virt_range(hfres):
          # Following Helgaker / Jørgensen / Olsen, p. 746:
          e_abij = orben[a] + orben[b] - orben[i] - orben[j]
          ga_aibj = eri[a, i, b, j] - eri[a, j, b, i]

          # Element ijab of the t2 amplitude
          t2[i,j,a,b] = - ga_aibj / e_abij

          # MP2 energy contribution:
          emp2 += 0.25 * (t2[i,j,a,b] * np.conj(ga_aibj)).real
  return emp2, t2

