#!/usr/bin/env python3

import molsturm
import water

params = {
  "atom_numbers": water.atom_numbers,
  "coords":       water.coords,
  #
  "basis_type":   "gaussian/libint",
  "basis_set":    "sto-3g",
}

print("Running calculation ...",end=" ",flush=True)
res = molsturm.hartree_fock(**params)
print("done")

print("\nFinal energies:")
molsturm.print_energies(res, indention=6*" ")

print("\nOrbital occupation:")
molsturm.print_mo_occupation(res,indention=6*" ")

molsturm.print_quote(res)
