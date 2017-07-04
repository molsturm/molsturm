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

  //! Restricted calculation or not
  bool restricted;

  /** \name SCF results */
  //@{
  //! Number of iterations the SCF took
  unsigned int n_iter = 0;

  //! Number of matrix applies needed
  unsigned int n_mtx_applies = 0;

  /** Frobenius norm of the final Pulay error matrix
   *  in the last step (which lead to convergence) */
  double final_error_norm = 0;

  /** Final change in total energy */
  double final_tot_energy_change = 0;

  /** Final change in one electron energy */
  double final_1e_energy_change = 0;
  //@}

  /** \name Energies */
  //@{
  //! Total SCF energy
  double energy_ground_state;

  //! Nuclear repulsion energy
  double energy_nuclear_repulsion;

  //! Electron-Nuclei attraction energy
  double energy_nuclear_attraction;

  //! Coulomb part of the  electron-electron interaction energy
  double energy_coulomb;

  //! Exchange part of the  electron-electron interaction energy
  double energy_exchange;

  //! Kinetic energy of the electron
  double energy_kinetic;
  //@}

  /** The expectation value of the total spin squared */
  double spin_squared;

  // Orbital coefficients
  std::vector<double> orbcoeff_bf;

  // Orbital energy
  std::vector<double> orben_f;

  // Fock matrix
  std::vector<double> fock_ff;

  // Overlap matrix
  std::vector<double> overlap_ff;

  // Matrix of all one electron integrals accumulated
  std::vector<double> hcore_ff;

  // Electron Repulsion Interaction tensor
  std::vector<double> eri_ffff;
};

}  // namespace iface
}  // namespace molsturm
