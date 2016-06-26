#pragma once
#include "IntegralOperatorBase.hh"
#include "RestrictedClosedIntegralOperator.hh"
#include <gscf/ScfStateBase.hh>
#include <linalgwrap/view.hh>

namespace molsturm {

template <typename IntegralOperator>
class ScfErrorLibrary {
public:
  static_assert(
        IsIntegralOperator<IntegralOperator>::value,
        "IntegralOperator needs to be derived off IntegralOperatorBase");

  typedef IntegralOperator operator_type;
  typedef typename operator_type::scalar_type scalar_type;
  typedef typename operator_type::stored_matrix_type matrix_type;
  typedef typename operator_type::size_type size_type;

  /** \brief Calculate the pulay error matrix in the space spanned by the
   *         basis functions.
   *
   * Constructs the lazy expression S * P * F - F * P * S, where
   * S is the overlap matrix from the state, F is the Restricted operator of the
   * state and P is the density computed from the alpha and beta values in F
   * and the eigenvalues of the state.
   */
  static auto pulay_error(const matrix_type& overlap_bb,
                          const matrix_type& coefficients_bf,
                          const operator_type& fock_bb)
        -> matrix_type {
    using namespace linalgwrap;

    static_assert(
          std::is_same<IntegralOperator,
                       RestrictedClosedIntegralOperator<matrix_type>>::value,
          "Currently this implementation only works for the "
          "RestrictedClosedIntegralOperator.");
    const size_type n_alpha = fock_bb.n_alpha();
    const size_type n_beta  = fock_bb.n_beta();

    // Occupied coefficients
    assert_size(n_alpha, n_beta);
    auto ca_bo = view::columns(coefficients_bf, range(n_alpha));

    // Density:
    auto pa_bb = ca_bo * view::transpose(ca_bo);

    // Return error expression, not evaluated
    // == S * P * F - F * P * S
    return static_cast<matrix_type>(view::view(overlap_bb) * pa_bb * fock_bb -
				    fock_bb * pa_bb * view::view(overlap_bb));
  }
};

}  // namespace molsturm
