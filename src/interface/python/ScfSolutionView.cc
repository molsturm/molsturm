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

#include "ScfSolutionView.hh"
#include "ExcInvalidParameters.hh"
#include <krims/ExceptionSystem.hh>

namespace molsturm {
namespace iface {

ScfSolutionView::ScfSolutionView(double* orben_sf_, int n_spin_e, int n_fock_e,
                                 double* orbcoeff_sbf_, int n_spin_c, int n_bas_c,
                                 int n_fock_c)
      : orbcoeff_sbf(orbcoeff_sbf_),
        orben_sf(orben_sf_),
        n_fock(static_cast<size_t>(n_fock_c)),
        n_bas(static_cast<size_t>(n_bas_c)),
        n_spin(static_cast<size_t>(n_spin_c)) {
  assert_throw(
        n_fock_e == n_fock_c,
        ExcInvalidParameters("Size mismatch in the number of orbitals in orben_f (" +
                             std::to_string(n_fock_e) + ") and orbcoeff_bf (" +
                             std::to_string(n_fock_c) + ")."));
  assert_throw(n_spin_e == n_spin_c,
               ExcInvalidParameters("Mismatch in number of spin components in orben_f (" +
                                    std::to_string(n_spin_e) + ") and orbcoeff_bf (" +
                                    std::to_string(n_spin_c) + ")."));
}

/** Obtain a copy as a lazyten::Eigensolution object */
lazyten::Eigensolution<scalar_type, vector_type> ScfSolutionView::as_eigensolution()
      const {

  // If we are unrestricted (2 spin components) then we need to build the block-diagonal
  // structure inside the eigensolution, which we pass onto the guess_external method
  // and hence to the initial Fock matrix. We do this by shifting the beta block by
  // n_fock and n_bas and placing it at the shifted location.
  lazyten::Eigensolution<double, lazyten::SmallVector<double>> esolution;
  auto& evalues  = esolution.evalues();
  auto& evectors = esolution.evectors();

  evalues.resize(n_spin * n_fock);
  std::copy(orben_sf, orben_sf + n_spin * n_fock, evalues.begin());
  evectors = lazyten::MultiVector<lazyten::SmallVector<double>>(n_spin * n_bas,
                                                                n_spin * n_fock);

  for (size_t s = 0; s < n_spin; ++s) {
    for (size_t b = 0; b < n_bas; ++b) {
      const size_t sb = s * n_bas + b;
      for (size_t f = 0; f < n_fock; ++f) {
        evectors[s * n_fock + f][s * n_bas + b] = orbcoeff_sbf[sb * n_fock + f];
      }  // f
    }    // b
  }      // s

  return esolution;
}

/** Obtain a copy of the coefficients of a particular spin component as a MultiVector */
lazyten::MultiVector<vector_type> ScfSolutionView::orbcoeff_bf(size_t spin) const {
  assert_greater(spin, n_spin);

  lazyten::MultiVector<vector_type> block(n_bas, n_fock);
  for (size_t b = 0; b < n_bas; ++b) {
    const size_t sb = spin * n_bas + b;
    for (size_t f = 0; f < n_fock; ++f) {
      block[f][b] = orbcoeff_sbf[sb * n_fock + f];
    }  // f
  }    // b

  return block;
}

/** Set values from a lazyten::Eigensolution object */
void ScfSolutionView::set_from_eigensolution(
      const lazyten::Eigensolution<scalar_type, vector_type>& eigensolution) {
  assert_throw(eigensolution.evalues().size() == n_spin * n_fock,
               krims::ExcSizeMismatch(eigensolution.evalues().size(), n_spin * n_fock));
  assert_throw(
        eigensolution.evectors().n_vectors() == n_spin * n_fock,
        krims::ExcSizeMismatch(eigensolution.evectors().n_vectors(), n_spin * n_fock));
  assert_throw(eigensolution.evectors().n_elem() == n_spin * n_bas,
               krims::ExcSizeMismatch(eigensolution.evectors().n_elem(), n_spin * n_bas));
  auto& evalues  = eigensolution.evalues();
  auto& evectors = eigensolution.evectors();

  std::copy(evalues.begin(), evalues.end(), orben_sf);
  for (size_t s = 0; s < n_spin; ++s) {
    for (size_t b = 0; b < n_bas; ++b) {
      const size_t sb = s * n_bas + b;
      for (size_t f = 0; f < n_fock; ++f) {
        orbcoeff_sbf[sb * n_fock + f] = evectors[s * n_fock + f][s * n_bas + b];
      }  // f
    }    // b
  }      // s

#ifdef DEBUG
  // Check the other blocks are zero.
  for (size_t s1 = 0; s1 < n_spin; ++s1) {
    for (size_t s2 = 0; s2 < n_spin; ++s2) {
      if (s1 == s2) continue;

      for (size_t b = 0; b < n_bas; ++b) {
        for (size_t f = 0; f < n_fock; ++f) {
          assert_internal(evectors[s1 * n_fock + f][s2 * n_bas + b] == 0);
        }  // f
      }    // b

    }   // s2
  }     // s1
#endif  // DEBUG
}

}  // namespace iface
}  // namespace molsturm
