#pragma once
#include <gscf/PulayDiisScfKeys.hh>

namespace molsturm {

struct IopDiisScfKeys : public gscf::PulayDiisScfKeys {
  /** Maximal total energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_tot_energy_change;

  /** Maximal 1e energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_1e_energy_change;
};
}  // namespace molsturm
