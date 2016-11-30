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
   * state and P is the density computed from the alpha and beta values in F
   * and the eigenvalues of the state.
   */
  template <typename OverlapMatrix>
  static auto pulay_error(
        const OverlapMatrix& overlap_bb,
        const linalgwrap::MultiVector<const vector_type>& coefficients_bf,
        const operator_type& fock_bb) -> matrix_type {
    const size_type n_alpha = fock_bb.n_alpha();
#ifdef DEBUG
    const size_type n_beta = fock_bb.n_beta();
#endif
    assert_dbg(n_alpha == n_beta, krims::ExcNotImplemented());
    assert_dbg(overlap_bb.is_symmetric(), linalgwrap::ExcMatrixNotSymmetric());
    assert_dbg(fock_bb.is_symmetric(), linalgwrap::ExcMatrixNotSymmetric());

    // Occupied coefficients
    auto ca_bo = coefficients_bf.subview({0, n_alpha});

    // Form first products (Factor 2 since alpha == beta)
    //    -- O(2*n_bas*n_bas*n_occ)
    auto Sca_bo = 2. * overlap_bb * ca_bo;
    auto Fca_bo = fock_bb * ca_bo;

    // Form the antisymmetric outer product sum and return it.
    // The idea is
    // S * P * F - F * P * S == S * C * C^T * F - F * C * C^T * S
    //                       == (S*C) * (F*C)^T - (F*C) * (S*C)^T
    //  -- O(n_bas*n_bas*n_occ)
    return outer_prod_sum(Sca_bo, Fca_bo) - outer_prod_sum(Fca_bo, Sca_bo);
  }
};
}  // namespace molsturm
