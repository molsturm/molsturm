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
#include <gint/sturmian/atomic/NlmBasis.hh>
#include <iostream>
#include <molsturm/MolecularSystem.hh>
#include <sstream>
#include <string>

namespace hf {
using gint::sturmian::atomic::NlmBasis;

struct args_type {
  // The system we model
  molsturm::MolecularSystem system;
  bool restricted = true;

  // Basis
  std::string basis_type;

  // Sturmians
  bool sturmian = false;
  double k_exp  = 1.0;
  int n_max     = 0;
  int l_max     = 0;
  int m_max     = 0;
  NlmBasis nlm_basis;

  // Gaussians
  bool gaussian         = false;
  std::string basis_set = "<not avail>";

  // Convergence
  size_t max_iter           = 25;
  double error              = 5e-7;
  size_t diis_size          = 4;
  size_t n_eigenpairs       = 0;
  std::string eigensolver   = "auto";
  std::string guess_esolver = "auto";
  std::string guess_method  = "hcore";
};

/** Write the content of args_type to a stream */
std::ostream& operator<<(std::ostream& o, const args_type& args);

/** Quick and dirty function to parse a string to a different type.
 *  Return false if not possible */
template <typename T>
bool str_to_type(const std::string& in, T& out) {
  return static_cast<bool>(std::stringstream(in) >> out);
}

/** \brief Quick and dirty function to parse the commandline arguments
 *
 * \returns true if all is fine, else false
 */
bool parse_args(int argc, char** argv, args_type& parsed);

}  // namespace hf
