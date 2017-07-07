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

/** Parameters which are available from the python interface.
 *  All parameters prefixed with internal_ are internal and will not
 *  be exposed to the user or available for the user
 */
struct Parameters {
  static const int all;

  //
  // System and setup
  //
  // The system we model
  int charge          = 0;
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

  bool internal_restricted_set_by_user = false;
  /** Use restricted or unrestricted fock operator
   * Note: This parameter is only used if
   * internal_restricted_set_by_user is also set to true.
   * */
  bool restricted = false;

  // Basis
  std::string basis_type = "";

  // Sturmians
  double k_exp = 0.0;
  int n_max    = 0;
  int l_max    = all;
  int m_max    = all;
  // TODO Phase out n_max, l_max and m_max and only transfer
  //      the nlm_basis array at all to the c++ side from python
  //      See also molsturm/scf_guess/_impl.py

  /** Transferred as a linearised array, i.e.
   *  (n1,l1,m1,n2,l2,m2, ... )
   */
  std::vector<int> nlm_basis = {};

  // Gaussians
  std::string basis_set = "";

  //
  // SCF and convergence
  //
  size_t max_iter           = 25;
  double conv_tol           = 5e-7;
  size_t diis_size          = 4;
  size_t n_eigenpairs       = 10000;
  std::string eigensolver   = "auto";
  std::string guess_esolver = "auto";
  std::string guess         = "hcore";

  // Printing
  bool print_iterations = false;

  // Special guess parameters for guess = external
  std::vector<double> guess_external_orben_f{};
  std::vector<double> guess_external_orbcoeff_bf{};

  //
  // Influence on export back to python
  //
  // Compute the repulsion integrals as a full tensor in the molecular orbital basis.
  // The underlying ao2mo transformation is rather slow.
  bool export_repulsion_integrals = false;

  // Export the fock matrix
  bool export_fock_matrix = false;

  // Export the matrix of all one electron terms combined
  bool export_hcore_matrix = true;

  // Export the overlap matrix (in MO basis)
  bool export_overlap_matrix = false;
};

}  // namespace iface
}  // namespace molsturm
