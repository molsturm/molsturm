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

import molsturm

params = {
  "atom_numbers": [4],
  "coords":       [[0,0,0]],
  #
  "basis_type":   "atomic/cs_reference_pc",
  "k_exp":        1.3,
  #
  "eigensolver":   "lapack",
  "guess_esolver": "lapack",
}

def run(**extra):
  params.update(extra)
  res = molsturm.hartree_fock(**params)
  molsturm.print_convergence_summary(res)
  molsturm.print_energies(res)
  molsturm.print_mo_occupation(res)
  molsturm.print_quote(res)
  return res

if __name__ == "__main__":
  run()
