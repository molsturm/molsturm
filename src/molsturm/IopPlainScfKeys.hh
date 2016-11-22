#pragma once
#include <gscf/PlainScfKeys.hh>

namespace molsturm {

struct IopPlainScfKeys : public gscf::PlainScfKeys {
  /** Maximal frobenius norm of the Pulay error matrix for convergence
   * (Type: real_type) */
  static const std::string max_error_norm;

  /** Maximal total energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_tot_energy_change;

  /** Maximal 1e energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_1e_energy_change;

  /** Should individual orbital energies be printed each iteration
   * (Type: bool */
  static const std::string print_orbital_energies;
};
}  // namespace molsturm
