//
// Copyright (C) 2017 by the molsturm authors
//
// This file is part of molsturm.
//
// molsturm is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// molsturm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with molsturm. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include <gscf/ScfBaseKeys.hh>

namespace molsturm {

struct IopScfKeys : public gscf::ScfBaseKeys {
  /** Pulay error norm below which DIIS acceleration is switched on */
  static const std::string diis_startup_error_norm;

  /** Iteration count after which DIIS acceleration is switched on
   *  (Regardless of the error norm */
  static const std::string diis_startup_iter;

  /** Pulay error norm below which DIIS is again switched off */
  static const std::string diis_shutdown_error_norm;

  /** Maximal total energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_tot_energy_change;

  /** Maximal 1e energy change between two cycles for convergence
   * (Type: real_type) */
  static const std::string max_1e_energy_change;

  /** Print the progress during the solve
   *  (Type: bool)
   *
   *  Will add ScfMsgType::IterationProcess to the verbosity parameter.
   */
  static const std::string print_iterations;

  /** Verbosity for the scf solver
   * (Type: ScfMsgType) */
  static const std::string verbosity;
};
}  // namespace molsturm
