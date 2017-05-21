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

molsturm.print_convergence_summary(res)
molsturm.print_energies(res)
molsturm.print_mo_occupation(res)
molsturm.print_quote(res)
