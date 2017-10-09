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

#ifdef SWIG
// clang-format off

// Setting the molecular structure
%apply (long* IN_ARRAY1, int DIM1) {(long* atom_numbers, int n_atoms_an)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2)
                {(double* coords, int n_atoms_c, int three_c)};

// Setting the nlm vectors for coulomb-sturmians
%apply (long* IN_ARRAY2, int DIM1, int DIM2) {(long* nlm, int n_nlm, int three_n)};

// clang-format on
#else
#include <gint/Structure.hh>
#include <krims/GenMap.hh>
#endif  // SWIG

namespace molsturm {
namespace iface {

#ifndef SWIG
/** The keys used by the system submap of the ScfParameters */
struct SystemKeys {
  /** The structure to use. Type: gint::Structure */
  static const std::string structure;

  /** The number of alpha electrons. Type: size_t */
  static const std::string n_alpha;

  /** The number of beta electrons. Type: size_t */
  static const std::string n_beta;
};
#endif  // SWIG

/** ScfParameters which are available from the python interface.
 *
 * Contains the following submaps:
 *    - guess     Will be used as the guess parameters i.e. to setup the scf guess
 *    - scf       Will be used as the scf parameters to be passed to the actual SCF solver
 *    - integrals Controls the integral backend, i.e. the basis function types which are
 *                   used and the parameters for these.
 *    - system    The chemical system to be modelled.
 *                   Contains the description of the nuclei and electrons.
 */
struct ScfParameters : public krims::GenMap {
  //
  // Getters and setters for the GenMap
  //
  /** Set a parameter in one of the GenMap objects */
  template <typename T>
  void update(const std::string& key, T value) {
    krims::GenMap::update(key, value);
  }

  //
  // Special parameter setters
  //
  /** Construct and set an gint::sturmian::NlmBasis object. */
  void update_nlm_basis(const std::string& key, long* nlm, int n_nlm, int three_n);

  /** Interpret this type string as a gint::OrbitalType enum and set the value */
  void update_orbital_type(const std::string& key, std::string type);

  /** Construct and set a gint::Structure object. */
  void update_structure(const std::string& key, long* atom_numbers, int n_atoms_an,
                        double* coords, int n_atoms_c, int three_c);
};

#ifdef SWIG
  // clang-format off
%extend ScfParameters {
  %template(update_bool)    update<bool>;
  %template(update_int)     update<int>;
  %template(update_size_t)  update<size_t>;
  %template(update_scalar)  update<double>;
  %template(update_string)  update<std::string>;
} // %extend
// clang-format on
#endif

}  // namespace iface
}  // namespace molsturm
