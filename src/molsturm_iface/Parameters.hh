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
#include <string>
#include <vector>

namespace molsturm {
namespace iface {

struct Parameters {
  static const int all;

  // The system we model
  int charge = 0;
  size_t multiplicity = 0;

  /** The coordinates of the atoms, transferred as
   *  a linearised array, i.e. the elements have the
   *  meaning
   *  [x1,y1,z1,x2,y2,z2, ... ]
   *
   *  Atomic units are assumed.
   */
  std::vector<double> coords = {};

  /** The atomic numbers of the atoms */
  std::vector<int> atom_numbers = {};

  /** The atom labels */
  std::vector<std::string> atoms = {};

  /** Use restricted or unrestricted fock operator
   * Note: This parameter is only used if restricted_set
   * is also set to true.
   * */
  bool restricted = false;
  bool restricted_set_by_user = false;

  // Basis
  std::string basis_type = "";

  // Sturmians
  double k_exp = 1.0;
  int n_max = 0;
  int l_max = all;
  int m_max = all;

  /** Transferred as a linearised array, i.e.
   *  (n1,l1,m1,n2,l2,m2, ... )
   */
  std::vector<int> nlm_basis = {};

  // Gaussians
  std::string basis_set = "";

  // Convergence
  size_t max_iter = 25;
  double error = 5e-7;
  size_t diis_size = 4;
  size_t n_eigenpairs = 10000;
  std::string eigensolver = "auto";
  std::string guess_esolver = "auto";
  std::string guess_method = "hcore";

  // Printing
  bool print_iterations = false;

  //
  // Influence on what is computed
  //
  // Compute the repulsion integrals as a full tensor in the molecular orbital basis.
  // The underlying ao2mo transformation is rather slow.
  bool enable_repulsion_integrals = true;
};

}  // namespace iface
}  // namespace molsturm
