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

#include "ScfParameters.hh"
#include "ExcInvalidParameters.hh"
#include "gint/config.hh"
#include <gint/IntegralLookupKeys.hh>
#include <gint/OrbitalType.hh>
#include <gint/sturmian/atomic/NlmBasis.hh>
#include <lazyten/Base/Solvers/Eigensolution.hh>
#include <lazyten/SmallMatrix.hh>
#include <limits>
#include <molsturm/scf_guess.hh>

namespace molsturm {
namespace iface {

const std::string SystemKeys::structure = "structure";
const std::string SystemKeys::n_alpha   = "n_alpha";
const std::string SystemKeys::n_beta    = "n_beta";

void ScfParameters::update_nlm_basis(const std::string& key, long* nlm, int n_nlm,
                                     int three_n) {
  using gint::sturmian::atomic::Nlm;
  using gint::sturmian::atomic::NlmBasis;

  assert_throw(three_n == 3,
               ExcInvalidParameters("The nlm array needs to be (n_bas, 3), i.e. the "
                                    "column dimension needs to be 3."));
  const int qnum_max = 1000;

  const size_t n_nlm_size = static_cast<size_t>(n_nlm);
  NlmBasis nlm_basis;
  nlm_basis.reserve(n_nlm_size);
  for (size_t i = 0; i < n_nlm_size; ++i) {
    Nlm tple{static_cast<int>(nlm[3 * i]), static_cast<int>(nlm[3 * i + 1]),
             static_cast<int>(nlm[3 * i + 2])};

    // Check values (mostly because we narrow long -> int)
    assert_throw(tple.n > 0 && tple.n < qnum_max,
                 ExcInvalidParameters("Values for n need to be in the range [1, " +
                                      std::to_string(qnum_max - 1) + "] and not " +
                                      std::to_string(tple.n) + "."));
    assert_throw(
          tple.l >= 0 && tple.l < tple.n,
          ExcInvalidParameters("Values for l need to be in the range [0, n-1] == [0, " +
                               std::to_string(tple.n - 1) + "] and not " +
                               std::to_string(tple.l) + "."));
    assert_throw(
          std::abs(tple.m) <= tple.l,
          ExcInvalidParameters("Values for m need to be in the range [-l, l] == [" +
                               std::to_string(-tple.l) + ", " + std::to_string(tple.l) +
                               "] and not " + std::to_string(tple.m) + "."));

    nlm_basis.push_back(std::move(tple));
  }

  krims::GenMap::update(key, std::move(nlm_basis));
}

/** Interpret this type string as a gint::OrbitalType enum and set the value */
void ScfParameters::update_orbital_type(const std::string& key, std::string type) {
  using gint::OrbitalType;

  const OrbitalType otype = [&type] {
    if (type == std::string("real_molecular")) {
      return OrbitalType::REAL_MOLECULAR;
    } else if (type == std::string("real_atomic")) {
      return OrbitalType::REAL_ATOMIC;
    } else if (type == std::string("complex_molecular")) {
      return OrbitalType::COMPLEX_MOLECULAR;
    } else if (type == std::string("complex_atomic")) {
      return OrbitalType::COMPLEX_ATOMIC;
    } else {
      assert_throw(false, ExcInvalidParameters(
                                "Invalid value for set_integral_param_orbital_type"));
      return OrbitalType::REAL_MOLECULAR;
    }
  }();

  krims::GenMap::update(key, otype);
}

/** Construct and set a gint::Structure object. */
void ScfParameters::update_structure(const std::string& key, long* atom_numbers,
                                     int n_atoms_an, double* coords, int n_atoms_c,
                                     int three_c) {
  size_t n_atoms = static_cast<size_t>(n_atoms_an);
  assert_throw(n_atoms_an == n_atoms_c,
               krims::ExcSizeMismatch(n_atoms, static_cast<size_t>(n_atoms_c)));
  assert_throw(three_c == 3, krims::ExcSizeMismatch(static_cast<size_t>(three_c), 3ul));
  static_assert(std::is_same<gint::real_type, double>::value,
                "gint::real_type is expected to be double right now.");

  // Build the structure object from the plain arrays:
  gint::Structure structure;
  structure.reserve(n_atoms);
  for (size_t i = 0; i < n_atoms; ++i) {
    std::array<double, 3> arr_coords{
          {coords[3 * i + 0], coords[3 * i + 1], coords[3 * i + 2]}};
    assert_throw(atom_numbers[i] > 0, krims::ExcTooLarge<long>(0, atom_numbers[i]));
    structure.emplace_back(atom_numbers[i], std::move(arr_coords));
  }

  krims::GenMap::update(key, std::move(structure));
}

}  // namespace iface
}  // namespace molsturm
