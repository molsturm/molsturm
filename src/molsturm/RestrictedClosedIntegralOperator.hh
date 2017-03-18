#pragma once
#include "IntegralOperatorBase.hh"
#include "IntegralTermContainer.hh"
#include "MolecularSystem.hh"
#include "energy_nuclear_repulsion.hh"
#include <functional>
#include <linalgwrap/Constants.hh>
#include <linalgwrap/LazyMatrixExpression.hh>
#include <linalgwrap/LazyMatrix_i.hh>
#include <map>

namespace molsturm {

// Open-shell: Have a simple specialisation of the eigensolver for now
//             and use adl (check that adl gets used properly)
//

// TODO Have an UHF version and a Restricted OHF version as well.
/** Class representing an restricted closed-shell integral operator */
template <typename StoredMatrix>
class RestrictedClosedIntegralOperator final : public IntegralOperatorBase<StoredMatrix> {
 public:
  typedef IntegralOperatorBase<StoredMatrix> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  /** Type of the coefficients and the pointer to the coefficients
   *  used inside the operator */
  typedef const linalgwrap::MultiVector<vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! Type of a integral term
  typedef typename IntegralTermContainer<stored_matrix_type>::int_term_type int_term_type;

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
   *        If the provided basis set behind the integrals has
   *        \t nbas basis functions, then this matrix should have
   *        the dimensionality \t nbas x \t nfock, where \t nfock is the
   *        number of fock operator eigenstates we calculate.
   * \param nalpha           Number of alpha electrons
   * \param nbeta            Number of beta electrons
   *                         (has to be equal to nalpha)
   */
  RestrictedClosedIntegralOperator(IntegralTermContainer<StoredMatrix> integral_terms,
                                   const MolecularSystem& system);

  /** Return the number of rows of the matrix */
  size_type n_rows() const override { return m_operator.n_rows(); }

  /** Return the number of columns of the matrix */
  size_type n_cols() const override { return m_operator.n_cols(); }

  /** Return an element of the matrix */
  scalar_type operator()(size_type row, size_type col) const override;

  bool has_transpose_operation_mode() const override {
    return m_operator.has_transpose_operation_mode();
  }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   *
   *  Loosely speaking we perform
   *  \[ M = c_M \cdot M + (A^{mode})_{rowrange,colrange} \]
   *  where
   *    - rowrange = [start_row, start_row+in.n_rows() ) and
   *    - colrange = [start_col, start_col+in.n_cols() )
   *
   * More details can be found in the same function in
   * LazyMatrixExpression
   */
  void extract_block(
        stored_matrix_type& M, const size_type start_row, const size_type start_col,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_M = linalgwrap::Constants<scalar_type>::zero) const override {
    m_operator.extract_block(M, start_row, start_col, mode, c_this, c_M);
  }

  /** \brief Compute the Matrix-Multivector application -- generic version
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   *
   * \note Whenever the virtual apply method is overwritten, this method
   * should be implemented as well as it assures that conversion to
   * MultiVector<MutableMemoryVector_i<scalar_type>> can actually occur
   * automatically.
   */
  template <typename VectorIn, typename VectorOut,
            linalgwrap::mat_vec_apply_enabled_t<RestrictedClosedIntegralOperator,
                                                VectorIn, VectorOut>...>
  void apply(const linalgwrap::MultiVector<VectorIn>& x,
             linalgwrap::MultiVector<VectorOut>& y,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_y = linalgwrap::Constants<scalar_type>::zero) const {
    m_operator.apply(x, y, mode, c_this, c_y);
  }

  /** \brief Compute the Matrix-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  void apply(
        const linalgwrap::MultiVector<
              const linalgwrap::MutableMemoryVector_i<scalar_type>>& x,
        linalgwrap::MultiVector<linalgwrap::MutableMemoryVector_i<scalar_type>>& y,
        const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
        const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
        const scalar_type c_y = linalgwrap::Constants<scalar_type>::zero) const override {
    m_operator.apply(x, y, mode, c_this, c_y);
  }

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   *
   * See LazyMatrixExpression for more details
   */
  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const linalgwrap::Transposed mode = linalgwrap::Transposed::None,
             const scalar_type c_this = linalgwrap::Constants<scalar_type>::one,
             const scalar_type c_out =
                   linalgwrap::Constants<scalar_type>::zero) const override {
    m_operator.mmult(in, out, mode, c_this, c_out);
  }

  /** Clone the matrix */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(
          new RestrictedClosedIntegralOperator<StoredMatrix>(*this));
  }

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   *
   * It expects the new coefficients under the parameter key
   * returned by scf_update_key()
   */
  void update(const krims::GenMap& map) override;

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   */
  void update(const coefficients_ptr_type& coefficients);

  /** Obtain the coefficients which were used to update the state of this operator,
   *  i.e. this returns the coefficients which have been most recently passed
   *  to this operator via update()
   */
  const coefficients_type& coefficients() const { return *m_coefficients_ptr; }

  /** Return the update key for the solver */
  const std::string& scf_update_key() const override { return m_update_key; }

  /* \name Access to energies of the individual electronic terms */
  ///@{
  /** Return the map from the ids of the integral terms specified on
   *  construction time to the energy value of that corresponding
   *  term.
   *
   * \note The nuclear repulsion energy can be obtained using the
   *       function call energy_nuclear_repulsion()
   */
  const std::map<gint::IntegralIdentifier, scalar_type>& energies() const {
    return m_energies;
  }

  /** Return the sum of the energies of all one-electron terms */
  scalar_type energy_1e_terms() const override;

  /** Return the sum of the energies of all two-electron terms */
  scalar_type energy_2e_terms() const override {
    return m_energies.at(m_coul.id()) + m_energies.at(m_exchge.id());
  }

  /** Return the nuclear repulsion component */
  scalar_type energy_nuclear_repulsion() const { return m_nuclear_repulsion_energy; }

  /** Return the total energy
   *
   * This means the *electronic* and the *nuclear* component. */
  scalar_type energy_total() const {
    return energy_1e_terms() + energy_2e_terms() + energy_nuclear_repulsion();
  }

  /** Return the 1e terms */
  std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
  terms_1e() const;

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of alpha
   * spin. */
  std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
  terms_alpha() const;

  // TODO it would be really great to have a better interface
  // for the function above and the one below.
  // E.g. one could use the block-diagonal matrix guys

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of beta
   * spin. */
  std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
  terms_beta() const {
    return terms_alpha();
  }
  ///@}

  krims::Range<size_t> indices_subspace(gscf::OrbitalSpace osp) const override {
    using gscf::OrbitalSpace;
    switch (osp) {
      case OrbitalSpace::OCC_ALPHA:
      /* intentional fall-through */
      case OrbitalSpace::OCC_BETA:
        return {0, m_n_alpha};
      case OrbitalSpace::VIRT_ALPHA:
      /* intentional fall-through */
      case OrbitalSpace::VIRT_BETA:
        return {m_n_alpha, n_rows()};
    }
  }

 private:
  /** Update the state of the lazy matrix operator terms and rebuild
   *  the operator m_operator from them.
   */
  void update_operator(const coefficients_ptr_type& coeff_bf_ptr);

  /** Recompute the energies from the current state of the terms
   * (First call update_operator on an update! )
   **/
  void update_energies(const coefficients_ptr_type& coeff_bf_ptr);

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
   * This operator is rebuild on each update call with the updated
   * matrices and thus is available for performing operations on
   * all terms together
   * */
  linalgwrap::LazyMatrixSum<stored_matrix_type> m_operator;

  /** Number of alpha electrons */
  const size_type m_n_alpha;

  /** The current electronic energies */
  std::map<gint::IntegralIdentifier, scalar_type> m_energies;

  /** The nuclear repulsion energy */
  scalar_type m_nuclear_repulsion_energy;

  /** The current coefficients the fock operator is initialised to */
  coefficients_ptr_type m_coefficients_ptr;

  //! Key used for updating the state.
  const std::string m_update_key = "evec_coefficients";
};

//
// --------------------------------------------------------------------
//

template <typename StoredMatrix>
RestrictedClosedIntegralOperator<StoredMatrix>::RestrictedClosedIntegralOperator(
      IntegralTermContainer<StoredMatrix> integral_terms, const MolecularSystem& system)
      : m_coeff_1e{std::move(integral_terms.coefficients_1e)},
        m_terms_1e{std::move(integral_terms.integral_terms_1e)},
        m_coeff_coul{integral_terms.coefficient_coulomb},
        m_coul{std::move(integral_terms.coulomb_term)},
        m_coeff_exchge{integral_terms.coefficient_exchange},
        m_exchge{std::move(integral_terms.exchange_term)},
        m_operator{},
        m_n_alpha{system.n_alpha},
        m_energies{} {
  using namespace linalgwrap;

  // Check that alpha is equal to beta
  assert_equal(system.n_beta, system.n_alpha);

  // Check that number of terms and number of coefficients agrees:
  assert_size(m_terms_1e.size(), m_coeff_1e.size());

  // Check operator size and zero energy terms
  const size_type op_size = m_coul.n_rows();

  auto itterm = std::begin(m_terms_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm) {
    assert_size(itterm->n_cols(), op_size);
    assert_size(itterm->n_rows(), op_size);

    m_energies.insert(std::make_pair(itterm->id(), Constants<scalar_type>::zero));
  }

  assert_size(m_coul.n_cols(), op_size);
  assert_size(m_coul.n_rows(), op_size);
  m_energies.insert(std::make_pair(m_coul.id(), Constants<scalar_type>::zero));

  assert_size(m_exchge.n_cols(), op_size);
  assert_size(m_exchge.n_rows(), op_size);
  m_energies.insert(std::make_pair(m_exchge.id(), Constants<scalar_type>::zero));

  m_nuclear_repulsion_energy = molsturm::energy_nuclear_repulsion(system.structure);

  // Initialise as a core hamiltonian (only 1e terms)
  const coefficients_ptr_type zero_coefficients =
        std::make_shared<coefficients_type>(op_size, m_n_alpha);
  update_operator(zero_coefficients);
}

template <typename StoredMatrix>
typename RestrictedClosedIntegralOperator<StoredMatrix>::scalar_type
RestrictedClosedIntegralOperator<StoredMatrix>::operator()(size_type row,
                                                           size_type col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  return m_operator(row, col);
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update(
      const coefficients_ptr_type& coeff_bf_ptr) {
  // Only update if we got new coefficients:
  if (coeff_bf_ptr == m_coefficients_ptr) return;
  m_coefficients_ptr = coeff_bf_ptr;

  update_operator(coeff_bf_ptr);
  update_energies(coeff_bf_ptr);
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update(const krims::GenMap& map) {
  if (map.exists(m_update_key)) {
    auto coeff_bf_ptr =
          static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(m_update_key));
    update(coeff_bf_ptr);
  }
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update_operator(
      const coefficients_ptr_type& coeff_bf_ptr) {
  // Here we assume restricted closed-shell HF, since we have just one
  // type of coefficients (alpha) to do the update.
  //
  // In fact we hence have twice the density represented by these coefficients
  // for the Coulomb terms. This is dealt with in the construction of the
  // operator below, where we provide an extra factor of 2 for the Coulomb term.

  // Build multivector of occupied alpha orbitals
  auto ca_bo_ptr = std::make_shared<const linalgwrap::MultiVector<const vector_type>>(
        coeff_bf_ptr->subview({0, m_n_alpha}));
  krims::GenMap occa_map{{int_term_type::update_key_coefficients, ca_bo_ptr}};

  // Initialise an empty operator
  m_operator = linalgwrap::LazyMatrixSum<StoredMatrix>();

  // Set up one-electron terms of m_fock:
  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    // Update the term, then add it to the operator:
    itterm->update(occa_map);
    m_operator += (*itcoeff) * (*itterm);
  }

  // We are closed shell --> just add it with an extra factor of 2
  // since alpha and beta density are identical and both
  // contribute to the coulomb operator
  m_coul.update(occa_map);
  m_operator += 2. * m_coeff_coul * m_coul;

  // Add exchange term (-1 should be done in coefficients):
  m_exchge.update(occa_map);
  m_operator += m_coeff_exchge * m_exchge;
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update_energies(
      const coefficients_ptr_type& coeff_bf_ptr) {
  // Here we assume restricted closed-shell HF, since we have just one
  // type of coefficients (alpha) to do the update.
  //
  // In fact we hence have twice the density represented by these coefficients
  // for the Coulomb terms. This is dealt with in the construction of the
  // operator above, where we provide an extra factor of 2 for the Coulomb term.

  // Build multivector of occupied alpha orbitals
  auto ca_bo_ptr = std::make_shared<const linalgwrap::MultiVector<const vector_type>>(
        coeff_bf_ptr->subview({0, m_n_alpha}));
  const auto& ca_bo = *ca_bo_ptr;
  krims::GenMap occa_map{{int_term_type::update_key_coefficients, ca_bo_ptr}};

  // Calculate the energies of the 1e terms:
  //
  // We use :
  // Energy = \tr (  \sum_{k \in occ} (C^{k})^T Op C^{k} )
  //        = \tr (  \sum_{k \in occ} (Op C^{k}) (C^{k})^T
  //        = \tr (  \sum_{k \in occ} outer_product(Op C^{k}, C^{k})
  //
  // where Op is an operator whose energy we want to compute and
  // C^{k} is the coefficient of the kth orbital
  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    // Calculate energy. Factor 2 because of alpha == beta
    const scalar_type energy = 2. * trace(outer_prod_sum((*itterm) * ca_bo, ca_bo));
    // TODO rethink this once we have the stored multivectors

    // Scale energy appropriately and set it in map:
    m_energies[itterm->id()] = (*itcoeff) * energy;
  }

  // Calculate the energies of the 2e terms:
  {
    // Coulomb energy:  Let J[P] be the Coulomb matrix build from the density P
    //     energy =   tr( P^{total} J[P^{total}] )
    //            =   tr( ( P^\alpha + P^\beta ) J[P^\alpha + P^beta] )
    //            = 2*tr( P^\alpha  2*J[P^\alpha] )
    //            = 4*tr( P^\alpha J[P^\alpha] )
    const scalar_type energy_coul = 4. * trace(outer_prod_sum(m_coul * ca_bo, ca_bo));

    // Exchange energy:  Let K[P] be the Exchange matrix build from the density P
    //     energy =   tr( P^\alpha K[P^\alpha] ) + tr( P^\beta K[P^\beta] )
    //            = 2 tr (P^\alpha K[P^\alpha] )
    const scalar_type energy_exchge = 2. * trace(outer_prod_sum(m_exchge * ca_bo, ca_bo));

    // Scale energy appropriately and set it in map:
    // 0.5 because those are 2e term and we need to avoid double counting.
    m_energies[m_coul.id()] = 0.5 * m_coeff_coul * energy_coul;
    m_energies[m_exchge.id()] = 0.5 * m_coeff_exchge * energy_exchge;
  }
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
std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
RestrictedClosedIntegralOperator<StoredMatrix>::terms_1e() const {
  using linalgwrap::LazyMatrixProduct;
  std::map<gint::IntegralIdentifier, LazyMatrixProduct<StoredMatrix>> ret;

  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    ret.insert(std::make_pair(itterm->id(), (*itcoeff) * (*itterm)));
  }

  return ret;
}

template <typename StoredMatrix>
std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
RestrictedClosedIntegralOperator<StoredMatrix>::terms_alpha() const {
  auto ret = terms_1e();
  ret.insert(std::make_pair(m_coul.id(), m_coeff_coul * m_coul));
  ret.insert(std::make_pair(m_exchge.id(), m_coeff_exchge * m_exchge));
  return ret;
}

}  // namespace molsturm
