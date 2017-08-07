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
#include "IntegralOperatorBase.hh"
#include "RestrictedClosedIntegralOperator.hh"
#include <gscf/ScfStateBase.hh>

namespace molsturm {

template <typename IntegralOperator>
class ScfErrorLibrary {
 public:
  static_assert(IsIntegralOperator<IntegralOperator>::value,
                "IntegralOperator needs to be derived off IntegralOperatorBase");

  typedef IntegralOperator operator_type;
  typedef typename operator_type::scalar_type scalar_type;
  typedef typename operator_type::stored_matrix_type matrix_type;
  typedef typename matrix_type::vector_type vector_type;
  typedef typename operator_type::size_type size_type;

  /** \brief Calculate the pulay error matrix in the space spanned by the
   *         basis functions.
   *
   * Constructs the expression S * P * F - F * P * S, where
   * S is the overlap matrix from the state, F is the Restricted operator of the
   * state and P is the density computed from the coefficients, alpha and beta
   * values in F.
   */
  template <typename OverlapMatrix>
  static matrix_type pulay_error(const operator_type& fock_bb,
                                 const OverlapMatrix& overlap_bb,
                                 bool accumulate_both = true) {
    auto occ_a = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
    auto occ_b = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA);

#ifdef DEBUG
    typedef typename OverlapMatrix::real_type real_type;
    const real_type tol = 100 * lazyten::Constants<real_type>::default_tolerance;
    assert_dbg(overlap_bb.is_symmetric(tol), lazyten::ExcMatrixNotSymmetric());
    assert_dbg(fock_bb.is_symmetric(tol), lazyten::ExcMatrixNotSymmetric());
#endif

    // Lambda to compute the error for a block (alpha-alpha or beta-beta)
    auto compute_error = [&fock_bb, &overlap_bb](krims::Range<size_t> occ) {
      auto cocc_bo = fock_bb.coefficients().subview(occ);

      // Form first products
      //    -- O(2*n_bas*n_bas*n_occ)
      auto Sc_bo = overlap_bb * cocc_bo;
      auto Fc_bo = fock_bb * cocc_bo;

      // Form the antisymmetric outer product sum
      //
      // The idea is
      // S * P * F - F * P * S == S * C * C^T * F - F * C * C^T * S
      //                       == (S*C) * (F*C)^T - (F*C) * (S*C)^T
      //  -- O(n_bas*n_bas*n_occ)
      return outer_prod_sum(Sc_bo, Fc_bo) - outer_prod_sum(Fc_bo, Sc_bo);
    };

    // Closed-shell case
    if (occ_a == occ_b) return 2. * compute_error(occ_a);

    // Compute individual errors in blocks
    auto error_alpha = compute_error(occ_a);
    auto error_beta  = compute_error(occ_b);

    if (accumulate_both) {
      // Sum the errors of the beta electrons onto the
      // appropriate alpha electrons (leaving out the
      // parts where we have less betas then alphas)
      for (size_t i = 0; i < error_beta.n_rows(); ++i) {
        for (size_t j = 0; j < error_beta.n_cols(); ++j) {
          error_alpha(i, j) += error_beta(i, j);
        }
      }
      return error_alpha;
    }

    // Build a block-diagonal matrix
    matrix_type error(error_alpha.n_rows() + error_beta.n_rows(),
                      error_beta.n_cols() + error_alpha.n_cols());
    for (size_t i = 0; i < error_alpha.n_rows(); ++i) {
      for (size_t j = 0; j < error_alpha.n_cols(); ++j) {
        error(i, j) += error_alpha(i, j);
      }
    }

    for (size_t i = 0; i < error_beta.n_rows(); ++i) {
      for (size_t j = 0; j < error_beta.n_cols(); ++j) {
        const size_t ib = error_alpha.n_rows() + i;
        const size_t jb = error_alpha.n_cols() + j;
        error(ib, jb) += error_beta(i, j);
      }
    }
    return error;
  }
};
}  // namespace molsturm
