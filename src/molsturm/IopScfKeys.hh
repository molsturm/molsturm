#pragma once
#include <gscf/ScfBaseKeys.hh>

namespace molsturm {

struct IopScfKeys : public gscf::ScfBaseKeys {
  /** Maximal total energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_tot_energy_change;

  /** Maximal 1e energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_1e_energy_change;

  /** Verbosity for the scf solver
   * (Type: ScfMsgType) */
  static const std::string verbosity;
};
}  // namespace molsturm
