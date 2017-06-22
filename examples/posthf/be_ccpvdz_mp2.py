#!/usr/bin/env python3
import molsturm
import molsturm.posthf

params = {
  "atom_numbers": [4],
  "coords":       [[0,0,0]],
  #
  "basis_type":   "gaussian/libint",
  "basis_set":    "cc-pvdz",
  #
  "print_iterations": True,
  #
  "error":          1e-10,
  #
  "export_repulsion_integrals": True,
  "export_fock_matrix": True,
}

res = molsturm.hartree_fock(**params)
molsturm.print_convergence_summary(res)
molsturm.print_energies(res)
molsturm.print_mo_occupation(res)

res_adc = molsturm.posthf.mp2(res)
print("MP2 energy", res_adc["energy_mp2"])
print("tot energy", res_adc["energy_ground_state"])

molsturm.print_quote(res)
