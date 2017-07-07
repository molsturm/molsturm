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
#include "IntegralTermContainer.hh"
#include <functional>
#include <linalgwrap/Constants.hh>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SubscriptionPointer.hh>
#include <linalgwrap/view.hh>
#include <map>

namespace molsturm {

// This was a try for a generalisation of the RestrictedClosedIntegralOperator
// to a general restricted operator.
template <typename StoredMatrix>
class RestrictedIntegralOperator : public linalgwrap::LazyMatrix_i<StoredMatrix> {
 public:
  typedef linalgwrap::LazyMatrix_i<StoredMatrix> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  //! Type of a integral term
  typedef typename IntegralTermContainer<stored_matrix_type>::int_term_type int_term_type;

  // TODO this is really bad, since it is a type which is only defined in
  // the detail namespace
  //! Type of a scale view to an integral term
  typedef linalgwrap::view::detail::ScaleView<int_term_type> scaled_int_term_type;

  /** \name Construct a Fock/Kohn-Sham operator for a restricted calculation
   *
   * It takes an IntegralTerms object, which contains information about all
   * the integral terms and the coefficients between them.
   *
   * I.e. if one passes it a list of integrals representing the projected
   * one-electron terms and the list of integrals representing the projected
   * two-electron terms with the coefficient 1 between them, it will represent
   * the default Fock operator for restricted Hartree-Fock.
   *
   * \param integral_terms   The IntegralTerms object.
   * \param initial_guess_bf The initial guess matrix.
   *        If the provided basis set behind the integrals has
   *        \t nbas basis functions, then this matrix should have
   *        the dimensionality \t nbas x \t nfock, where \t nfock is the
   *        number of fock operator eigenstates we calculate.
   * \param nalpha           Number of alpha electrons
   * \param nbeta            Number of beta electrons
   *                         (has to be less or equal to nalpha)
   */
  RestrictedIntegralOperator(IntegralTermContainer<StoredMatrix> integral_terms,
                             const stored_matrix_type& initial_guess_bf,
                             size_type n_alpha, size_type n_beta);

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
  void update(const linalgwrap::GenMap& map) override;

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

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of alpha
   * spin. */
  std::map<std::string, scaled_int_term_type> terms_alpha() const;

  // TODO it would be really great to have a better interface
  // for the function above and the one below.
  // E.g. one could use the block-diagonal matrix guys

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of beta
   * spin. */
  std::map<std::string, scaled_int_term_type> terms_beta() const;
  ///@}

  /** Get the number of alpha electrons */
  size_type n_alpha() const;

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

  // TODO
  // Think about a way to have the alpha and beta exchange terms
  // and alpha and beta coulomb terms
  // in parallel such that the overhead is as small as possible
  // and can potentially even be resolved at compile time.
  //
  // (e.g. templated inner helper-class, which takes a template
  // argument, which states whether alpha == beta)

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

  /** Number of alpha electrons. We assume alpha <= beta */
  const size_type m_n_beta;

  /** The current energies */
  std::map<std::string, scalar_type> m_energies;
};

//
// --------------------------------------------------------------------
//

template <typename StoredMatrix>
RestrictedIntegralOperator<StoredMatrix>::RestrictedIntegralOperator(
      IntegralTermContainer<StoredMatrix> integral_terms,
      const stored_matrix_type& initial_guess_bf, size_type n_alpha, size_type n_beta)
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
  assert_greater_equal(n_beta, n_alpha);

#ifdef DEBUG
  // Operator size:
  size_type op_size = integral_terms.coulomb_term.n_rows();
#endif

  // Check that number of terms and number of coefficients agrees:
  assert_size(m_terms_1e.size(), m_coeff_1e.size());

  // Set up one-electron terms of m_fock:
  auto itterm  = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    // Are matrices square and of the expected size?
    assert_dbg(itterm->n_rows() == itterm->n_cols(), ExcMatrixNotSquare());
    assert_dbg(itterm->n_rows() == op_size, ExcMatrixNotSquare());

    // Create a view to the term, scale it and add it to the fock:
    m_operator += view::scale(*itterm, *itcoeff);

    // Set all energies to zero:
    m_energies.insert(itterm->id(), Constants<scalar_type>::zero);
  }

  // Add coulomb term:
  if (n_alpha == n_beta) {
    // Closed shell --> just add it with an extra factor of 2
    // since alpha and beta density are identical and both
    // contribute to the coulomb operator
    m_operator += view::scale(m_coul, 2. * m_coeff_coul);
  } else {
    // TODO Restricted open shell --> not implemented.
    assert_implemented(false);
  }
  m_energies.insert(m_coul.id(), Constants<scalar_type>::zero);

  // Add exchange term (-1 should be done in coefficients):
  if (n_alpha == n_beta) {
    m_operator += view::scale(m_exchge, m_coeff_exchge);
  } else {
    // TODO Do not know the RO-HF case atm
    assert_implemented(false);
  }
  m_energies.insert(m_exchge.id(), Constants<scalar_type>::zero);

  // TODO note for later:
  // For the uhf we have to be more careful here,
  // since all terms but the K matrices are identical, but those
  // are in fact different. So we best be very careful when designing
  // the block diagonal matrices and adjusting the terms here.
  // or do everything implicitly in this operator structure.

  // Build the first proper fock matrix and its energies:
  update_state_and_energies(initial_guess_bf);
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::size_type
RestrictedIntegralOperator<StoredMatrix>::n_rows() const {
  return m_operator.n_rows();
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::size_type
RestrictedIntegralOperator<StoredMatrix>::n_cols() const {
  return m_operator.n_cols();
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::scalar_type
RestrictedIntegralOperator<StoredMatrix>::operator()(size_type row, size_type col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  return m_operator(row, col);
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::stored_matrix_type
      RestrictedIntegralOperator<StoredMatrix>::operator*(
            const stored_matrix_type& in) const {
  assert_size(m_operator.n_cols(), in.n_rows());
  return m_operator * in;
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::lazy_matrix_expression_ptr_type
RestrictedIntegralOperator<StoredMatrix>::clone() const {
  return lazy_matrix_expression_ptr_type(
        new RestrictedIntegralOperator<StoredMatrix>(*this));
}

template <typename StoredMatrix>
void RestrictedIntegralOperator<StoredMatrix>::update(const linalgwrap::GenMap& map) {
  // The coefficient key we look for:
  const std::string coeff_key("evec_coefficients");

  if (map.exists(coeff_key)) {
    // We have new coefficients
    const auto& coeff = map.at<stored_matrix_type>(coeff_key);
    update_state_and_energies(coeff);
  }
}

template <typename StoredMatrix>
void RestrictedIntegralOperator<StoredMatrix>::update_state_and_energies(
      const stored_matrix_type& coefficients_bf) {
  using namespace linalgwrap;

  //
  // Occupied alpha and beta states in parameter maps:
  //
  auto ca_bo = view::columns(coefficients_bf, range(m_n_alpha));
  auto cb_bo = view::columns(coefficients_bf, range(m_n_beta));

  const std::string occ_coeff_key("coefficients_occupied");

  GenMap occa_map;
  occa_map.update(occ_coeff_key,
                  make_subscription(ca_bo, "GenMap in RestrictedOperator"));

  GenMap occb_map;
  occb_map.update(occ_coeff_key,
                  make_subscription(cb_bo, "GenMap in RestrictedOperator"));

  //
  // Update 1e part, J and K == Update the operator.
  //
  // TODO Do not know the RO-HF case atm
  assert_implemented(m_n_beta == m_n_alpha);

  // Here we assume resticted closed-shell HF, since we have just one
  // type of coefficients (alpha) to do the update.
  //
  // In fact we hence have twice the density represented by these coefficients
  // for the Coulomb terms. This is dealt with in the constructor, where we
  // provide an extra factor of 2 for the Coulomb term.
  m_operator.update(occa_map);

  //
  // Calculate the energies of the 1e terms:
  //
  auto itterm  = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    // TODO Do not know the RO-HF case atm
    assert_implemented(m_n_beta == m_n_alpha);

    // Calculate energy. Factor 2 because of alpha == beta
    scalar_type energy = 2. * (view::transpose(ca_bo) * (*itterm) * ca_bo).trace();

    // Scale energy appropriately and set it in map:
    m_energies[itterm->id()] = (*itcoeff) * energy;
  }

  //
  // Calculate the energies of the 2e terms:
  //
  {
    // TODO Do not know the RO-HF case atm
    assert_implemented(m_n_beta == m_n_alpha);

    // Coulomb: Factor 2 because of alpha == beta
    scalar_type energy_coul = 2. * (view::transpose(ca_bo) * m_coul * ca_bo).trace();

    // Exchange: Factor 2 because of alpha == beta
    scalar_type energy_exchge = 2. * (view::transpose(ca_bo) * m_exchge * ca_bo).trace();

    // Scale energy appropriately and set it in map:
    // 0.5 because those are 2e term and we need to avoid double counting.
    m_energies[m_coul.id()]   = 0.5 * m_coeff_coul * energy_coul;
    m_energies[m_exchge.id()] = 0.5 * m_coeff_exchge * energy_exchge;
  }
}

template <typename StoredMatrix>
const std::map<std::string,
               typename RestrictedIntegralOperator<StoredMatrix>::scalar_type>&
RestrictedIntegralOperator<StoredMatrix>::energies() const {
  return m_energies;
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::scalar_type
RestrictedIntegralOperator<StoredMatrix>::energy_1e_terms() const {
  // m_energies contains the correct values by the ids
  //
  // We use the collection integral_terms_1e to obtain the ids of the 1e terms
  // than we sum them all

  scalar_type sum = linalgwrap::Constants<scalar_type>::zero;
  for (auto& term : m_terms_1e) {
    sum += m_energies[term.id()];
  }
  return sum;
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::scalar_type
RestrictedIntegralOperator<StoredMatrix>::energy_2e_terms() const {
  return m_energies[m_coul.id()] + m_energies[m_exchge.id()];
}

template <typename StoredMatrix>
std::map<std::string,
         typename RestrictedIntegralOperator<StoredMatrix>::scaled_int_term_type>
RestrictedIntegralOperator<StoredMatrix>::terms_alpha() const {
  using namespace linalgwrap;

  assert_implemented(m_n_beta == m_n_alpha);

  std::map<std::string, scaled_int_term_type> ret;

  // Insert 1e terms:
  auto itterm  = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    ret.insert(itterm->id(), view::scale(*itterm, *itcoeff));
  }

  // Insert 2e terms:
  ret.insert(m_coul.id(), view::scale(m_coul, m_coeff_coul));
  ret.insert(m_exchge.id(), view::scale(m_exchge, m_coeff_exchge));

  return ret;
}

template <typename StoredMatrix>
std::map<std::string,
         typename RestrictedIntegralOperator<StoredMatrix>::scaled_int_term_type>
RestrictedIntegralOperator<StoredMatrix>::terms_beta() const {
  assert_implemented(m_n_beta == m_n_alpha);
  return terms_alpha();
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::size_type
RestrictedIntegralOperator<StoredMatrix>::n_alpha() const {
  return m_n_alpha;
}

template <typename StoredMatrix>
typename RestrictedIntegralOperator<StoredMatrix>::size_type
RestrictedIntegralOperator<StoredMatrix>::n_beta() const {
  return m_n_beta;
}

}  // namespace molsturm
