#pragma once
#include <gint/Integral.hh>
#include <linalgwrap/Constants.hh>

namespace molsturm {

template <typename StoredMatrix>
struct IntegralTermContainer {
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;

  //! The type to use for the integral terms.
  typedef typename gint::Integral<stored_matrix_type> int_term_type;

  /** The one electron terms */
  std::vector<int_term_type> integral_terms_1e;

  /** Set of coefficients for the one electron integral terms */
  std::vector<scalar_type> coefficients_1e;

  /** The coulomb term to use */
  int_term_type coulomb_term;

  /** The coefficient for the coulomb term */
  scalar_type coefficient_coulomb;

  /** The exchange term to use */
  int_term_type exchange_term;

  /** The coefficient for the HF exchange term */
  scalar_type coefficient_exchange;

  /** \brief Setup the class with the given integral terms.
   *  \note By default all coefficients, but the exchange coefficient, are set
   * to one. The exchange term coefficient is set to -1.
   *  This is done in order to make sure that by default this data structure is
   *  ready for a Hartree-Fock calculation.
   */
  IntegralTermContainer(std::vector<int_term_type> integral_terms_1e_,
                        int_term_type coulomb_term_, int_term_type exchange_term_)
        : integral_terms_1e{std::move(integral_terms_1e_)},
          coefficients_1e{std::vector<scalar_type>(
                integral_terms_1e.size(), linalgwrap::Constants<scalar_type>::one)},
          coulomb_term{std::move(coulomb_term_)},
          coefficient_coulomb{linalgwrap::Constants<scalar_type>::one},
          exchange_term{std::move(exchange_term_)},
          coefficient_exchange{-linalgwrap::Constants<scalar_type>::one} {}
};

}  // namespace molsturm
