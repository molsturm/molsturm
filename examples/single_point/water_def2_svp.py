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

molsturm.print_convergence_summary(res)
molsturm.print_energies(res)
molsturm.print_mo_occupation(res)
molsturm.print_quote(res)
