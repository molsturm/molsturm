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
#include <gint/IntegralUpdateKeys.hh>

namespace molsturm {

// TODO Derive of lazymatrixsum ... it reduces the boilerplate code a lot.

/** Class representing an restricted closed-shell integral operator */
template <typename StoredMatrix>
class RestrictedClosedIntegralOperator : public IntegralOperatorBase<StoredMatrix>,
                                         public lazyten::LazyMatrix_i<StoredMatrix> {
 public:
  typedef IntegralOperatorBase<StoredMatrix> base_type;
  typedef typename base_type::coefficients_ptr_type coefficients_ptr_type;
  typedef typename base_type::coefficients_type coefficients_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename lazyten::LazyMatrix_i<StoredMatrix>::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  constexpr bool restricted() const { return true; }

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
   * for a restricted closed-shell Hartree-Fock calculation.
   *
   * \param integral_terms   The IntegralTerms object.
   *        If the provided basis set behind the integrals has
   *        \t nbas basis functions, then this matrix should have
   *        the dimensionality \t nbas x \t nfock, where \t nfock is the
   *        number of fock operator eigenstates we calculate.
   * \param system   Molecular system use for the calculation
   */
  RestrictedClosedIntegralOperator(IntegralTermContainer<StoredMatrix> integral_terms,
                                   const MolecularSystem& system)
        : base_type{integral_terms, system}, m_operator{} {
    // Check that alpha is equal to beta
    assert_equal(system.n_beta, system.n_alpha);

    // Use only zero coefficients
    const size_t n_bas          = base_type::m_coul_adens.n_rows();
    coefficients_ptr_type zeros = std::make_shared<coefficients_type>(
          n_bas, base_type::m_n_alpha + base_type::m_n_beta);
    base_type::m_coefficients_ptr = zeros;
    update_state(zeros, /* include_energies= */ false);
  }

  /** Return the number of rows of the matrix */
  size_t n_rows() const override { return m_operator.n_rows(); }

  /** Return the number of columns of the matrix */
  size_t n_cols() const override { return m_operator.n_cols(); }

  /** Return an element of the matrix */
  scalar_type operator()(size_t row, size_t col) const override {
    return m_operator(row, col);
  }

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
  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col,
                     const lazyten::Transposed mode = lazyten::Transposed::None,
                     const scalar_type c_this       = 1,
                     const scalar_type c_M          = 0) const override {
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
            lazyten::mat_vec_apply_enabled_t<RestrictedClosedIntegralOperator, VectorIn,
                                             VectorOut>...>
  void apply(const lazyten::MultiVector<VectorIn>& x, lazyten::MultiVector<VectorOut>& y,
             const lazyten::Transposed mode = lazyten::Transposed::None,
             const scalar_type c_this = 1, const scalar_type c_y = 0) const {
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
        const lazyten::MultiVector<const lazyten::MutableMemoryVector_i<scalar_type>>& x,
        lazyten::MultiVector<lazyten::MutableMemoryVector_i<scalar_type>>& y,
        const lazyten::Transposed mode = lazyten::Transposed::None,
        const scalar_type c_this = 1, const scalar_type c_y = 0) const override {
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
             const lazyten::Transposed mode = lazyten::Transposed::None,
             const scalar_type c_this = 1, const scalar_type c_out = 0) const override {
    m_operator.mmult(in, out, mode, c_this, c_out);
  }

  /** Clone the matrix */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(
          new RestrictedClosedIntegralOperator<StoredMatrix>(*this));
  }

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of alpha
   * spin. */
  std::map<gint::IntegralIdentifier, lazyten::LazyMatrixProduct<StoredMatrix>>
  terms_alpha() const override final;

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of beta
   * spin. */
  std::map<gint::IntegralIdentifier, lazyten::LazyMatrixProduct<StoredMatrix>>
  terms_beta() const override final {
    return terms_alpha();
  }
  ///@}

  krims::Range<size_t> indices_orbspace(gscf::OrbitalSpace osp) const override final;

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   *
   * It expects the new coefficients under the parameter key
   * returned by scf_update_key()
   */
  void update(const krims::GenMap& map) override { base_type::update(map); }
  using base_type::update;

 private:
  //@{
  /** Types used to pass occupied coefficients to the exchange and coulomb parts */
  typedef const lazyten::MultiVector<const vector_type> cocc_type;
  typedef std::shared_ptr<const lazyten::MultiVector<const vector_type>> cocc_ptr_type;
  //@}

  /** Instruction to update energies and operator. Called by update in the base class */
  void update_state(const coefficients_ptr_type& coeff_bf_ptr,
                    bool include_energies) override final;

  /** Update the state of the lazy matrix operator terms and rebuild
   *  the operator m_operator from them.
   */
  void update_operator(const cocc_ptr_type& ca_bo_ptr);

  /** Recompute the energies from the current state of the terms
   * (First call update_operator on an update! )
   **/
  void update_energies(const cocc_ptr_type& ca_bo_ptr);

  /** The actual Operator matrix object as a LazyMatrixSum
   *
   * This operator is rebuild on each update call with the updated
   * matrices and thus is available for performing operations on
   * all terms together
   * */
  lazyten::LazyMatrixSum<stored_matrix_type> m_operator;
};

//
// --------------------------------------------------------------------
//

template <typename StoredMatrix>
krims::Range<size_t> RestrictedClosedIntegralOperator<StoredMatrix>::indices_orbspace(
      gscf::OrbitalSpace osp) const {
  using gscf::OrbitalSpace;
  assert_internal(base_type::m_coefficients_ptr != nullptr);
  const size_t n_orbs = base_type::m_coefficients_ptr->n_vectors();

  switch (osp) {
    case OrbitalSpace::OCC_ALPHA:
    /* intentional fall-through */
    case OrbitalSpace::OCC_BETA:
      return {0, base_type::m_n_alpha};
    case OrbitalSpace::VIRT_ALPHA:
    /* intentional fall-through */
    case OrbitalSpace::VIRT_BETA:
      return {base_type::m_n_alpha, n_orbs};
    default:
      return {0, 0};
  }
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update_state(
      const coefficients_ptr_type& coeff_bf_ptr, bool include_energies) {
  auto ca_bo_ptr = std::make_shared<const lazyten::MultiVector<const vector_type>>(
        coeff_bf_ptr->subview({0, base_type::m_n_alpha}));
  update_operator(ca_bo_ptr);

  if (include_energies) {
    update_energies(ca_bo_ptr);
  }
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update_operator(
      const cocc_ptr_type& ca_bo_ptr) {
  // Here we assume restricted closed-shell HF, since we have just one
  // type of coefficients (alpha) to do the update.
  //
  // In fact we hence have twice the density represented by these coefficients
  // for the Coulomb terms. This is dealt with in the construction of the
  // operator below, where we provide an extra factor of 2 for the Coulomb term.
  krims::GenMap occa_map{{gint::IntegralUpdateKeys::coefficients_occupied, ca_bo_ptr}};

  // Empty the operator:
  m_operator = lazyten::LazyMatrixSum<StoredMatrix>{};

  // Set up one-electron terms of m_fock:
  auto itterm  = std::begin(base_type::m_terms_1e);
  auto itcoeff = std::begin(base_type::m_coeff_1e);
  for (; itterm != std::end(base_type::m_terms_1e); ++itterm, ++itcoeff) {
    // Note: No update for one-electron terms needed, since no dependence on density
    m_operator += (*itcoeff) * (*itterm);
  }

  // We are closed shell --> just add it with an extra factor of 2
  // since alpha and beta density are identical and both
  // contribute to the coulomb operator
  base_type::m_coul_adens.update(occa_map);
  m_operator += 2. * base_type::m_coeff_coul * base_type::m_coul_adens;

  // Add exchange term (-1 should be done in coefficients):
  base_type::m_exchge_adens.update(occa_map);
  m_operator += base_type::m_coeff_exchge * base_type::m_exchge_adens;
}

template <typename StoredMatrix>
std::map<gint::IntegralIdentifier, lazyten::LazyMatrixProduct<StoredMatrix>>
RestrictedClosedIntegralOperator<StoredMatrix>::terms_alpha() const {
  // Note the factor 2 in front of the coulomb term once again comes from the fact
  // that we are closed-shell and hence only ever use the alpha density, but the
  // coulomb interaction is to both alpha and beta electrons.

  auto ret = base_type::terms_1e();
  ret.insert(std::make_pair(base_type::m_coul_adens.id(),
                            2. * base_type::m_coeff_coul * base_type::m_coul_adens));
  ret.insert(std::make_pair(base_type::m_exchge_adens.id(),
                            base_type::m_coeff_exchge * base_type::m_exchge_adens));
  return ret;
}

template <typename StoredMatrix>
void RestrictedClosedIntegralOperator<StoredMatrix>::update_energies(
      const cocc_ptr_type& ca_bo_ptr) {
  auto& coeff_coul   = base_type::m_coeff_coul;
  auto& coul_adens   = base_type::m_coul_adens;  // Coulomb from alpha density
  auto& coeff_exchge = base_type::m_coeff_exchge;
  auto& exchge_adens = base_type::m_exchge_adens;  // Exchange from alpha density
  auto& ca_bo        = *ca_bo_ptr;

  // Here we assume restricted closed-shell HF, since we have just one
  // type of coefficients (alpha) to do the update.
  //
  // In fact we hence have twice the density represented by these coefficients
  // for the Coulomb terms. This is dealt with in the construction of the
  // operator above, where we provide an extra factor of 2 for the Coulomb term.

  // Calculate the energies of the 1e terms:
  //
  // We use :
  // Energy = \tr (  \sum_{k \in occ} (C^{k})^T Op C^{k} )
  //        = \tr (  \sum_{k \in occ} (Op C^{k}) (C^{k})^T
  //        = \tr (  \sum_{k \in occ} outer_product(Op C^{k}, C^{k})
  //
  // where Op is an operator whose energy we want to compute and
  // C^{k} is the coefficient of the kth orbital
  auto itterm  = std::begin(base_type::m_terms_1e);
  auto itcoeff = std::begin(base_type::m_coeff_1e);
  for (; itterm != std::end(base_type::m_terms_1e); ++itterm, ++itcoeff) {
    // Calculate energy. Factor 2 because of alpha == beta
    const scalar_type energy = 2. * trace(outer_prod_sum((*itterm) * ca_bo, ca_bo));
    // TODO rethink this once we have the stored multivectors

    // Scale energy appropriately and set it in map:
    base_type::m_energies[itterm->id()] = (*itcoeff) * energy;
  }

  // Calculate the energies of the 2e terms:
  {
    // Coulomb energy:  Let J[P] be the Coulomb matrix build from the density P
    // both J[P] and P are block-diagonal, hence
    //     energy =   tr( P^\alpha J[P^{total}] ) + tr( P^\beta J[P^{total}] )
    //            =   tr( ( P^\alpha + P^\beta ) J[P^\alpha + P^beta] )
    //            = 2*tr( P^\alpha  2*J[P^\alpha] )
    //            = 4*tr( P^\alpha J[P^\alpha] )
    const scalar_type energy_coul = 4. * trace(outer_prod_sum(coul_adens * ca_bo, ca_bo));

    // Exchange energy:  Let K[P] be the Exchange matrix build from the density P
    // Again K[P] and P are block-diagonal
    //     energy =   tr( P^\alpha K[P^\alpha] ) + tr( P^\beta K[P^\beta] )
    //            = 2 tr (P^\alpha K[P^\alpha] )
    const scalar_type energy_exchge =
          2. * trace(outer_prod_sum(exchge_adens * ca_bo, ca_bo));

    // Scale energy appropriately and set it in map:
    // 0.5 because those are 2e term and we need to avoid double counting.
    base_type::m_energies[coul_adens.id()]   = 0.5 * coeff_coul * energy_coul;
    base_type::m_energies[exchge_adens.id()] = 0.5 * coeff_exchge * energy_exchge;
  }
}

}  // namespace molsturm
