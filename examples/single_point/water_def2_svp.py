#!/usr/bin/env python3

import molsturm
import water

params = {
  "atom_numbers": water.atom_numbers,
  "coords":       water.coords,
  #
  "basis_type":   "gaussian/libint",
  "basis_set":    "def2-svp",
  #
  "print_iterations": True,
}
res = molsturm.hartree_fock(**params)

print("\nFinal energies:")
molsturm.print_energies(res, indention=6*" ")

print("\nOrbital occupation:")
molsturm.print_mo_occupation(res,indention=6*" ")

molsturm.print_quote(res)
