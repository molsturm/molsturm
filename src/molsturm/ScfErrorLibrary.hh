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
                                 const OverlapMatrix& overlap_bb) {
    auto occ_a = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
    auto occ_b = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA);
    assert_implemented(occ_a == occ_b);

#ifdef DEBUG
    typedef typename OverlapMatrix::real_type real_type;
    const real_type tol = 100 * linalgwrap::Constants<real_type>::default_tolerance;
    assert_dbg(overlap_bb.is_symmetric(tol), linalgwrap::ExcMatrixNotSymmetric());
    assert_dbg(fock_bb.is_symmetric(tol), linalgwrap::ExcMatrixNotSymmetric());
#endif

    // Occupied coefficients for alpha
    auto ca_bo = fock_bb.coefficients().subview(occ_a);

    // Form first products for alpha (Factor 2 since alpha == beta)
    //    -- O(2*n_bas*n_bas*n_occ)
    auto Sca_bo = overlap_bb * ca_bo;
    auto Fca_bo = 2 * fock_bb * ca_bo;

    // Form the antisymmetric outer product sum for alpha
    //
    // The idea is
    // S * P * F - F * P * S == S * C * C^T * F - F * C * C^T * S
    //                       == (S*C) * (F*C)^T - (F*C) * (S*C)^T
    //  -- O(n_bas*n_bas*n_occ)
    return outer_prod_sum(Sca_bo, Fca_bo) - outer_prod_sum(Fca_bo, Sca_bo);
  }
};
}  // namespace molsturm
