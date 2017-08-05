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

import itertools as it

def occ_range(hfres):
  # Chain of the range of occupied alphas and the range
  # of occupied betas
  return it.chain(range(hfres["n_alpha"]),
                  range(hfres["n_orbs_alpha"],
                        hfres["n_orbs_alpha"] + hfres["n_beta"]))

def virt_range(hfres):
  # Chain of the range of virtual alphas and the range
  # of virtual betas
  return it.chain(range(hfres["n_alpha"], hfres["n_orbs_alpha"]),
                  range(hfres["n_orbs_alpha"] + hfres["n_beta"], 
                        hfres["n_orbs_alpha"] + hfres["n_orbs_beta"]))


