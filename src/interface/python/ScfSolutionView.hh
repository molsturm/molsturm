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
#include "config.hh"
#include <lazyten/Base/Solvers/Eigensolution.hh>

#ifdef SWIG
// clang-format off
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* orben_sf, int n_spin_e, int n_fock_e)};
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3)
                {(double* orbcoeff_sbf, int n_spin_c, int n_bas_c, int n_fock_c)};
// clang-format on
#endif  // SWIG

namespace molsturm {
namespace iface {
/** Structure representing a view into a SCF solution on a very low level (i.e. raw
 * pointers)
 *
 * Usually we use this to set data to arrays on the numpy side, which represent an SCF
 * solution.
 */
struct ScfSolutionView {
  /**
   *
   * \param p  ScfParameters object
   * \param orben_sf    On call contains the orbital energy
   *                    values of a potential
   *                    external guess for the SCF procedure.
   *                    (selected via the ScfParameters)
   *                    On exit contains the actual orbital energies
   *                    obtained by the SCF procedure.
   *                    The first index is the spin index (0 == alpha,
   *                    1 == beta) and the second index is the
   *                    orbital index. For restricted calculations
   *                    only one spin block should be specified.
   *  \param orbcoeff_sbf   On call contains the orbital coefficients of a
   *                        potential external guess for the SCF procedure.
   *                        on exit contains the converged coefficients.
   *                        The first index is spin again and for restricted
   *                        we only expect one block.
   */
  ScfSolutionView(double* orben_sf, int n_spin_e, int n_fock_e, double* orbcoeff_sbf,
                  int n_spin_c, int n_bas_c, int n_fock_c);

#ifndef SWIG
  double* orbcoeff_sbf;
  double* orben_sf;

  size_t n_fock;
  size_t n_bas;
  size_t n_spin;

  /** Obtain a copy as a lazyten::Eigensolution object */
  lazyten::Eigensolution<scalar_type, vector_type> as_eigensolution() const;

  /** Obtain a copy of the coefficients of a particular spin component as a MultiVector */
  lazyten::MultiVector<vector_type> orbcoeff_bf(size_t spin) const;

  /** Set values from a lazyten::Eigensolution object */
  void set_from_eigensolution(
        const lazyten::Eigensolution<scalar_type, vector_type>& eigensolution);
#endif  // SWIG
};

}  // namespace iface
}  // namespace molsturm
