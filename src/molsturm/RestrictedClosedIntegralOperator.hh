#pragma once
#include "IntegralOperatorBase.hh"
#include "IntegralTermContainer.hh"
#include <functional>
#include <linalgwrap/Constants.hh>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SubscriptionPointer.hh>
#include <linalgwrap/view.hh>
#include <map>

namespace molsturm {

// TODO Have an UHF version and a Restricted OHF version as well.
/** Class representing an restricted closed-shell integral operator */
template <typename StoredMatrix>
class RestrictedClosedIntegralOperator
      : public IntegralOperatorBase<StoredMatrix> {
public:
  typedef IntegralOperatorBase<StoredMatrix> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  //! Type of a integral term
  typedef typename IntegralTermContainer<stored_matrix_type>::int_term_type
        int_term_type;

  // TODO this is really bad, since it is a type which is only defined in
  // the detail namespace
  //! Type of a scale view to an integral term
  typedef linalgwrap::view::detail::ScaleView<const int_term_type>
        scaled_const_int_term_type;

  /** \name Construct a Fock/Kohn-Sham operator for a restricted closed-shell
   *  calculation
   *
   * It takes an IntegralTerms object, which contains information about all
   * the integral terms and the coefficients between them.
   *
   * I.e. if one passes it a list of integrals representing the projected
   * one-electron terms and the coulomb integral with coefficient 1 and
   * the exchange integral with coefficient -1 (The default of
   * IntegralTermContainer) it will represent the default Fock operator
   * for a restricted closed-shell Hartree-Fock calculation.list of integrals
   * representing the projected.
   *
   * \param integral_terms   The IntegralTerms object.
   * \param initial_guess_bf The initial guess matrix.
   *        If the provided basis set behind the integrals has
   *        \t nbas basis functions, then this matrix should have
   *        the dimensionality \t nbas x \t nfock, where \t nfock is the
   *        number of fock operator eigenstates we calculate.
   * \param nalpha           Number of alpha electrons
   * \param nbeta            Number of beta electrons
   *                         (has to be equal to nalpha)
   */
  RestrictedClosedIntegralOperator(
        IntegralTermContainer<StoredMatrix> integral_terms,
        const stored_matrix_type& initial_guess_bf, size_type n_alpha,
        size_type n_beta);

  /** Return the number of rows of the matrix */
  size_type n_rows() const override;

  /** Return the number of columns of the matrix */
  size_type n_cols() const override;

  /** Return an element of the matrix */
  scalar_type operator()(size_type row, size_type col) const override;

  /** Apply matrix to another matrix */
  stored_matrix_type operator*(const stored_matrix_type& in) const override;

  /** Clone the matrix */
  lazy_matrix_expression_ptr_type clone() const override;

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   *
   * It expects the new coefficients under the parameter key
   * "evec_coefficients"
   */
  void update(const linalgwrap::ParameterMap& map) override;

  /* \name Access to energies and individual terms */
  ///@{
  /** Return the map from the ids of the integral terms specified on
   *  construction time to the energy value of that corresponding
   *  term.
   */
  const std::map<std::string, scalar_type>& energies() const;

  /** Return the sum of the energies of all one-electron terms */
  scalar_type energy_1e_terms() const;

  /** Return the sum of the energies of all two-electron terms */
  scalar_type energy_2e_terms() const;

  /** Return the total energy */
  scalar_type energy_total() const;

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of alpha
   * spin. */
  std::map<std::string, scaled_const_int_term_type> terms_alpha() const;

  // TODO it would be really great to have a better interface
  // for the function above and the one below.
  // E.g. one could use the block-diagonal matrix guys

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of beta
   * spin. */
  std::map<std::string, scaled_const_int_term_type> terms_beta() const;
  ///@}

  /** Get the number of alpha electrons */
  size_type n_alpha() const;

  // TODO being able to access the above and below here feels wrong.
  // but we need it for the error calculation.

  /** Get the number of beta electrons */
  size_type n_beta() const;

private:
  /** Update the state of m_operator and calculate the new energies */
  void update_state_and_energies(const stored_matrix_type& coefficients_bf);

  /** Set of coefficients for the one electron integral terms */
  std::vector<scalar_type> m_coeff_1e;

  /** The one electron terms */
  std::vector<int_term_type> m_terms_1e;

  /** The coefficient for the coulomb term */
  scalar_type m_coeff_coul;

  /** The coulomb term to use */
  int_term_type m_coul;

  /** The coefficient for the HF exchange term */
  scalar_type m_coeff_exchge;

  /** The exchange term to use */
  int_term_type m_exchge;

  /** The actual Operator matrix object as a LazyMatrixSum
   *
   * Note that this is a sum made out of views to the integral
   * terms above. This makes sure, that the terms in this object
   * and the terms in the m_integral_terms object share exactly
   * the same state, as copying views is just creating another
   * view to exactly the same object.
   * */
  linalgwrap::LazyMatrixSum<stored_matrix_type> m_operator;

  /** Number of alpha electrons */
  const size_type m_n_alpha;

  /** The current energies */
  std::map<std::string, scalar_type> m_energies;
};

//
// --------------------------------------------------------------------
//

template <typename StoredMatrix>
RestrictedClosedIntegralOperator<StoredMatrix>::
      RestrictedClosedIntegralOperator(
            IntegralTermContainer<StoredMatrix> integral_terms,
            const stored_matrix_type& initial_guess_bf, size_type n_alpha,
            size_type n_beta)
      : m_coeff_1e{std::move(integral_terms.coefficients_1e)},
        m_terms_1e{std::move(integral_terms.integral_terms_1e)},
        m_coeff_coul{integral_terms.coefficient_coulomb},
        m_coul{std::move(integral_terms.coulomb_term)},
        m_coeff_exchge{integral_terms.coefficient_exchange},
        m_exchge{std::move(integral_terms.exchange_term)},
        m_operator{},
        m_n_alpha{n_alpha},
        m_energies{} {
  using namespace linalgwrap;

  // Check that alpha is greater-equal to beta:
  assert_equal(n_beta, n_alpha);

#ifdef DEBUG
  // Operator size:
  size_type op_size = integral_terms.coulomb_term.n_rows();
#endif

  // Check that number of terms and number of coefficients agrees:
  assert_size(m_terms_1e.size(), m_coeff_1e.size());

  // Set up one-electron terms of m_fock:
  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    // Are matrices square and of the expected size?
    assert_dbg(itterm->n_rows() == itterm->n_cols(), ExcMatrixNotSquare());
    assert_dbg(itterm->n_rows() == op_size, ExcMatrixNotSquare());

    // Create a view to the term, scale it and add it to the fock:
    m_operator += view::scale(*itterm, *itcoeff);

    // Set all energies to zero:
    m_energies.insert(
          std::make_pair(itterm->id(), Constants<scalar_type>::zero));
  }

  // We are closed shell --> just add it with an extra factor of 2
  // since alpha and beta density are identical and both
  // contribute to the coulomb operator
  m_operator += view::scale(m_coul, 2. * m_coeff_coul);
  m_energies.insert(std::make_pair(m_coul.id(), Constants<scalar_type>::zero));

  // Add exchange term (-1 should be done in coefficients):
  m_operator += view::scale(m_exchge, m_coeff_exchge);
  m_energies.insert(
        std::make_pair(m_exchge.id(), Constants<scalar_type>::zero));

  // Build the first proper fock matrix and its energies:
  update_state_and_energies(initial_guess_bf);
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::size_type
RestrictedClosedIntegralOperator<StoredMatrix>::n_rows() const {
  return m_operator.n_rows();
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::size_type
RestrictedClosedIntegralOperator<StoredMatrix>::n_cols() const {
  return m_operator.n_cols();
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::scalar_type
RestrictedClosedIntegralOperator<StoredMatrix>::operator()(
      size_type row, size_type col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  return m_operator(row, col);
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::stored_matrix_type
      RestrictedClosedIntegralOperator<StoredMatrix>::operator*(
            const stored_matrix_type& in) const {
  assert_size(m_operator.n_cols(), in.n_rows());
  return m_operator * in;
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<
      StoredMatrix>::lazy_matrix_expression_ptr_type
RestrictedClosedIntegralOperator<StoredMatrix>::clone() const {
  return lazy_matrix_expression_ptr_type(
        new RestrictedClosedIntegralOperator<StoredMatrix>(*this));
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update(
      const linalgwrap::ParameterMap& map) {
  // The coefficient key we look for:
  const std::string coeff_key("evec_coefficients");

  if (map.exists(coeff_key)) {
    // We have new coefficients
    const auto& coeff = map.at<stored_matrix_type>(coeff_key);
    update_state_and_energies(coeff);
  }
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update_state_and_energies(
      const stored_matrix_type& coefficients_bf) {
  using namespace linalgwrap;

  //
  // Occupied alpha states in parameter maps:
  //
  auto ca_bo = view::columns(coefficients_bf, range(m_n_alpha));

  const std::string occ_coeff_key("coefficients_occupied");

  // Maybe pass occupation number / occupation matrix and
  // full coefficients (as pointer) instead?
  ParameterMap occa_map;
  occa_map.update_copy(occ_coeff_key,
                       std::move(static_cast<stored_matrix_type>(ca_bo)));

  //
  // Update 1e part, J and K == Update the operator.
  //

  // Here we assume restricted closed-shell HF, since we have just one
  // type of coefficients (alpha) to do the update.
  //
  // In fact we hence have twice the density represented by these coefficients
  // for the Coulomb terms. This is dealt with in the constructor, where we
  // provide an extra factor of 2 for the Coulomb term.
  m_operator.update(occa_map);

  //
  // Calculate the energies of the 1e terms:
  //
  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    // Calculate energy. Factor 2 because of alpha == beta
    scalar_type energy =
          2. * (view::transpose(ca_bo) * (*itterm) * ca_bo).trace();

    // Scale energy appropriately and set it in map:
    m_energies[itterm->id()] = (*itcoeff) * energy;
  }

  //
  // Calculate the energies of the 2e terms:
  //
  {
    // Coulomb: Factor 2 because of alpha == beta
    scalar_type energy_coul =
          2. * (view::transpose(ca_bo) * m_coul * ca_bo).trace();

    // Exchange: Factor 2 because of alpha == beta
    scalar_type energy_exchge =
          2. * (view::transpose(ca_bo) * m_exchge * ca_bo).trace();

    // Scale energy appropriately and set it in map:
    // 0.5 because those are 2e term and we need to avoid double counting.
    m_energies[m_coul.id()] = 0.5 * m_coeff_coul * energy_coul;
    m_energies[m_exchge.id()] = 0.5 * m_coeff_exchge * energy_exchge;
  }
}

template <typename StoredMatrix>
const std::map<std::string, typename RestrictedClosedIntegralOperator<
                                  StoredMatrix>::scalar_type>&
RestrictedClosedIntegralOperator<StoredMatrix>::energies() const {
  return m_energies;
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::scalar_type
RestrictedClosedIntegralOperator<StoredMatrix>::energy_1e_terms() const {
  // m_energies contains the correct values by the ids
  //
  // We use the collection integral_terms_1e to obtain the ids of the 1e terms
  // than we sum them all

  scalar_type sum = linalgwrap::Constants<scalar_type>::zero;
  for (auto& term : m_terms_1e) {
    sum += m_energies.at(term.id());
  }
  return sum;
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::scalar_type
RestrictedClosedIntegralOperator<StoredMatrix>::energy_2e_terms() const {
  return m_energies.at(m_coul.id()) + m_energies.at(m_exchge.id());
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::scalar_type
RestrictedClosedIntegralOperator<StoredMatrix>::energy_total() const {
  return energy_1e_terms() + energy_2e_terms();
}

template <typename StoredMatrix>
std::map<std::string, typename RestrictedClosedIntegralOperator<
                            StoredMatrix>::scaled_const_int_term_type>
RestrictedClosedIntegralOperator<StoredMatrix>::terms_alpha() const {
  using namespace linalgwrap;
  std::map<std::string, scaled_const_int_term_type> ret;

  // Insert 1e terms:
  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    ret.insert(std::make_pair(itterm->id(), view::scale(*itterm, *itcoeff)));
  }

  // Insert 2e terms:
  ret.insert(std::make_pair(m_coul.id(), view::scale(m_coul, m_coeff_coul)));
  ret.insert(
        std::make_pair(m_exchge.id(), view::scale(m_exchge, m_coeff_exchge)));

  return ret;
}

template <typename StoredMatrix>
std::map<std::string, typename RestrictedClosedIntegralOperator<
                            StoredMatrix>::scaled_const_int_term_type>
RestrictedClosedIntegralOperator<StoredMatrix>::terms_beta() const {
  return terms_alpha();
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::size_type
RestrictedClosedIntegralOperator<StoredMatrix>::n_alpha() const {
  return m_n_alpha;
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::size_type
RestrictedClosedIntegralOperator<StoredMatrix>::n_beta() const {
  return m_n_alpha;
}

}  // namespace molsturm
