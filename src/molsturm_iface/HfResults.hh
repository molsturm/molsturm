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
#include <vector>

namespace molsturm {
namespace iface {

struct HfResults {
  unsigned int n_alpha = 0;
  unsigned int n_beta = 0;
  unsigned int n_bas = 0;

  unsigned int n_orbs_alpha = 0;
  unsigned int n_orbs_beta = 0;

  //! SCF convergence threshold
  double threshold = 0;

  // TODO ideals
  //       - number of iterations
  //       - final error
  //       - more energies

  //! Total SCF energy
  double energy_total;

  //! Nuclear repulsion energy
  double energy_nuclear_repulsion;

  //! Restricted calculation or not
  bool restricted;

  std::vector<double> coeff_fb;
  std::vector<double> orbital_energies_f;
  std::vector<double> fock_ff;
  std::vector<double> repulsion_integrals_ffff;
};

}  // namespace iface
}  // namespace molsturm
