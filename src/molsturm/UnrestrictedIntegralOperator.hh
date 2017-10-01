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
#include <lazyten/BlockDiagonalMatrix.hh>

namespace molsturm {

template <typename StoredMatrix,
          typename BlockType = lazyten::LazyMatrixSum<StoredMatrix>>
class UnrestrictedIntegralOperator : public IntegralOperatorBase<StoredMatrix>,
                                     public lazyten::BlockDiagonalMatrix<BlockType, 2> {
 public:
  typedef IntegralOperatorBase<StoredMatrix> base_type;
  typedef typename base_type::coefficients_ptr_type coefficients_ptr_type;
  typedef typename base_type::coefficients_type coefficients_type;
  typedef StoredMatrix stored_matrix_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::int_term_type int_term_type;

  typedef typename lazyten::BlockDiagonalMatrix<
        BlockType, 2>::lazy_matrix_expression_ptr_type lazy_matrix_expression_ptr_type;

  constexpr bool restricted() const { return false; }

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
   * for a unrestricted Hartree-Fock calculation.
   *
   * \param integral_terms   The IntegralTerms object.
   *        If the provided basis set behind the integrals has
   *        \t nbas basis functions, then this matrix should have
   *        the dimensionality \t nbas x \t nfock, where \t nfock is the
   *        number of fock operator eigenstates we calculate.
   * \param system   Molecular system use for the calculation
   */
  UnrestrictedIntegralOperator(IntegralTermContainer<StoredMatrix> integral_terms,
                               const MolecularSystem& system);

  /** Clone the matrix */
  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(
          new UnrestrictedIntegralOperator<StoredMatrix>(*this));
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
  terms_beta() const override final;
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

  /** The beta coulomb term
   * The coulombic part resulting from the alpha electron density
   * \note This is not the part which sits in the alpha-alpha block !
   */
  int_term_type m_coul_bdens;

  /** The beta exchange term
   * The coulombic 2e part resulting from the alpha electron density
   * \note This is not the part which sits in the alpha-alpha block !
   */
  int_term_type m_exchge_bdens;

  /** Instruction to update energies and operator. Called by update in the base class */
  void update_state(const coefficients_ptr_type& coeff_bf_ptr,
                    bool include_energies) override final;

  /** Update the state of the lazy matrix operator terms and rebuild the operator for
   *  alpha and beta blocks
   *
   *  \param ca_bo    Alpha occupied coefficients
   *  \param cb_bo    Beta occupied coefficients
   */
  void update_operator(const cocc_ptr_type& ca_bo_ptr, const cocc_ptr_type& cb_bo_ptr);

  /** Recompute the energies from the current state of the terms
   * (First call update_operator on an update! )
   **/
  void update_energies(const cocc_ptr_type& ca_bo_ptr, const cocc_ptr_type& cb_bo_ptr);

  /** Check that the coefficients have zeros at the expected places */
  void assert_coefficient_structure(const coefficients_ptr_type& coeff_bf_ptr) const;
};

//
// --------------------------------------------------------------------
//

template <typename StoredMatrix, typename BlockType>
UnrestrictedIntegralOperator<StoredMatrix, BlockType>::UnrestrictedIntegralOperator(
      IntegralTermContainer<StoredMatrix> integral_terms, const MolecularSystem& system)
      : base_type{integral_terms, system},
        // We need to initialise the BlockDiagonalMatrix with something,
        // which will later resemble the actual operator, so we will use
        // the first 1e term as a dummy here
        lazyten::BlockDiagonalMatrix<BlockType, 2>{
              {{BlockType{integral_terms.integral_terms_1e[0]},
                BlockType{integral_terms.integral_terms_1e[0]}}}},
        m_coul_bdens{integral_terms.coulomb_term},
        m_exchge_bdens{integral_terms.exchange_term} {
  // Use only zero coefficients
  const size_t n_bas      = base_type::m_coul_adens.n_rows();
  const size_t n_max_elec = std::max(base_type::m_n_alpha, base_type::m_n_beta);
  coefficients_ptr_type zeros =
        std::make_shared<coefficients_type>(2 * n_bas, 2 * n_max_elec);
  base_type::m_coefficients_ptr = zeros;
  update_state(zeros, /* include_energies= */ false);
}

template <typename StoredMatrix, typename BlockType>
krims::Range<size_t>
UnrestrictedIntegralOperator<StoredMatrix, BlockType>::indices_orbspace(
      gscf::OrbitalSpace osp) const {
  using gscf::OrbitalSpace;
  assert_internal(base_type::m_coefficients_ptr != nullptr);
  const size_t n_orbs = base_type::m_coefficients_ptr->n_vectors();
  assert_internal(n_orbs % 2 == 0);
  const size_t n_orbs_alpha = n_orbs / 2;
  const size_t n_alpha      = base_type::m_n_alpha;
  const size_t n_beta       = base_type::m_n_beta;

  switch (osp) {
    case OrbitalSpace::OCC_ALPHA:
      return {0, n_alpha};
    case OrbitalSpace::VIRT_ALPHA:
      return {n_alpha, n_orbs_alpha};
    case OrbitalSpace::OCC_BETA:
      return {n_orbs_alpha, n_orbs_alpha + n_beta};
    case OrbitalSpace::VIRT_BETA:
      return {n_orbs_alpha + n_beta, n_orbs};
    default:
      return {0, 0};
  }
}

template <typename StoredMatrix, typename BlockType>
void UnrestrictedIntegralOperator<StoredMatrix, BlockType>::update_state(
      const coefficients_ptr_type& coeff_bf_ptr, bool include_energies) {
  assert_coefficient_structure(coeff_bf_ptr);

  auto& coeff_bf            = *coeff_bf_ptr;
  const size_t n_orbs_alpha = coeff_bf_ptr->n_vectors() / 2;
  const size_t n_bas        = m_coul_bdens.n_rows();
  const size_t n_alpha      = base_type::m_n_alpha;
  const size_t n_beta       = base_type::m_n_beta;

  assert_throw(coeff_bf_ptr->n_elem() == 2 * n_bas,
               krims::ExcSizeMismatch(coeff_bf_ptr->n_elem(), 2 * n_bas));
  assert_throw(n_alpha + n_beta <= coeff_bf_ptr->n_vectors(),
               krims::ExcTooLarge<size_t>(n_alpha + n_beta, coeff_bf_ptr->n_vectors()));

  // Since the eigenvectors from the Block-Diagonal eigensolver are padded with zeros
  // in the alpha-beta and beta-alpha blocks, we need to remove that padding and get
  // back the actual coefficient values in only alpha-alpha and beta-beta.
  //
  // TODO Right now the only way to do this is by copying

  lazyten::MultiVector<vector_type> reduced_coeff_bf(n_bas, n_alpha + n_beta, false);
  // Copy alpha-alpha block
  for (size_t v = 0; v < n_alpha; ++v) {
    std::copy(coeff_bf[v].begin(), coeff_bf[v].begin() + n_bas,
              reduced_coeff_bf[v].begin());
  }

  // Copy beta-beta block
  for (size_t v = 0; v < n_beta; ++v) {
    const size_t v_beta = v + n_orbs_alpha;
    std::copy(coeff_bf[v_beta].begin() + n_bas, coeff_bf[v_beta].end(),
              reduced_coeff_bf[v + n_alpha].begin());
  }

  // Make alpha and beta occupied views
  auto ca_bo_ptr = std::make_shared<cocc_type>(reduced_coeff_bf.subview({0, n_alpha}));
  auto cb_bo_ptr = std::make_shared<cocc_type>(
        reduced_coeff_bf.subview({n_alpha, n_alpha + n_beta}));

  update_operator(ca_bo_ptr, cb_bo_ptr);
  if (include_energies) {
    update_energies(ca_bo_ptr, cb_bo_ptr);
  }
}

template <typename StoredMatrix, typename BlockType>
void UnrestrictedIntegralOperator<StoredMatrix, BlockType>::update_operator(
      const cocc_ptr_type& ca_bo_ptr, const cocc_ptr_type& cb_bo_ptr) {
  auto& blocka       = this->diag_blocks()[0];  // alpha-alpha block
  auto& blockb       = this->diag_blocks()[1];  // beta-beta block
  auto& coeff_coul   = base_type::m_coeff_coul;
  auto& coul_adens   = base_type::m_coul_adens;  // Coulomb from alpha density
  auto& coeff_exchge = base_type::m_coeff_exchge;
  auto& exchge_adens = base_type::m_exchge_adens;  // Exchange from alpha density

  // Clear the current terms
  blocka = lazyten::LazyMatrixSum<StoredMatrix>{};
  blockb = lazyten::LazyMatrixSum<StoredMatrix>{};

  // Set up one-electron terms:
  auto itterm  = std::begin(base_type::m_terms_1e);
  auto itcoeff = std::begin(base_type::m_coeff_1e);
  for (; itterm != std::end(base_type::m_terms_1e); ++itterm, ++itcoeff) {
    // Note: No update for one-electron terms needed, since no dependence on density
    blocka += (*itcoeff) * (*itterm);
    blockb += (*itcoeff) * (*itterm);
  }

  krims::GenMap occa_map{{gint::IntegralUpdateKeys::coefficients_occupied, ca_bo_ptr}};
  krims::GenMap occb_map{{gint::IntegralUpdateKeys::coefficients_occupied, cb_bo_ptr}};

  // Build and add the coulomb terms
  coul_adens.update(occa_map);
  m_coul_bdens.update(occb_map);
  blocka += coeff_coul * coul_adens + coeff_coul * m_coul_bdens;
  blockb += coeff_coul * coul_adens + coeff_coul * m_coul_bdens;

  // Build and add the exchange terms
  // Note that there is no exchange interaction of the alpha electrons with the alpha
  // electrons
  // and vice versa due to spin symmetry reasons
  // This is why we need the exchange matrix made from the beta density in the block for
  // alpha!
  exchge_adens.update(occa_map);
  m_exchge_bdens.update(occb_map);
  blocka += coeff_exchge * exchge_adens;
  blockb += coeff_exchge * m_exchge_bdens;

  assert_internal(blocka.n_rows() == blocka.n_cols());
  assert_internal(blockb.n_rows() == blockb.n_cols());
  assert_internal(blockb.n_rows() == blocka.n_rows());
}

template <typename StoredMatrix, typename BlockType>
std::map<gint::IntegralIdentifier, lazyten::LazyMatrixProduct<StoredMatrix>>
UnrestrictedIntegralOperator<StoredMatrix, BlockType>::terms_alpha() const {
  auto ret = base_type::terms_1e();
  lazyten::LazyMatrixProduct<StoredMatrix> total_coul(base_type::m_coul_adens +
                                                      m_coul_bdens);
  ret.insert(std::make_pair(base_type::m_coul_adens.id(),
                            base_type::m_coeff_coul * total_coul));
  ret.insert(std::make_pair(base_type::m_exchge_adens.id(),
                            base_type::m_coeff_exchge * base_type::m_exchge_adens));
  return ret;
}

template <typename StoredMatrix, typename BlockType>
std::map<gint::IntegralIdentifier, lazyten::LazyMatrixProduct<StoredMatrix>>
UnrestrictedIntegralOperator<StoredMatrix, BlockType>::terms_beta() const {
  auto ret = base_type::terms_1e();
  lazyten::LazyMatrixProduct<StoredMatrix> total_coul(base_type::m_coul_adens +
                                                      m_coul_bdens);
  ret.insert(std::make_pair(base_type::m_coul_adens.id(),
                            base_type::m_coeff_coul * total_coul));
  ret.insert(std::make_pair(m_exchge_bdens.id(),  //
                            base_type::m_coeff_exchge * m_exchge_bdens));
  return ret;
}

template <typename StoredMatrix, typename BlockType>
void UnrestrictedIntegralOperator<StoredMatrix, BlockType>::update_energies(
      const cocc_ptr_type& ca_bo_ptr, const cocc_ptr_type& cb_bo_ptr) {
  auto& coeff_coul   = base_type::m_coeff_coul;
  auto& coul_adens   = base_type::m_coul_adens;  // Coulomb from alpha density
  auto& coeff_exchge = base_type::m_coeff_exchge;
  auto& exchge_adens = base_type::m_exchge_adens;  // Exchange from alpha density
  const auto& ca_bo  = *ca_bo_ptr;
  const auto& cb_bo  = *cb_bo_ptr;

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
    // Calculate energies:
    const scalar_type energy_alpha = trace(outer_prod_sum((*itterm) * ca_bo, ca_bo));
    const scalar_type energy_beta  = trace(outer_prod_sum((*itterm) * cb_bo, cb_bo));
    // TODO Rethink this once we have the stored multivectors

    // Scale energy appropriately and set it in map:
    base_type::m_energies[itterm->id()] = (*itcoeff) * (energy_alpha + energy_beta);
  }

  // Calculate the energies of the 2e terms:
  {
    // Coulomb energy:  Let J[P] be the Coulomb matrix build from the density P
    // both J[P] and P are block-diagonal, hence
    //     energy =   tr( P^\alpha J[P^{total}] ) + tr( P^\beta J[P^{total}] )
    //            =   tr( P^\alpha J[P^\alpha] ) + tr( P^\alpha J[P^\beta] )
    //              + tr( P^\beta  J[P^\alpha] ) + tr( P^\beta  J[P^\beta] )
    const scalar_type energy_coul = trace(outer_prod_sum(coul_adens * ca_bo, ca_bo)) +
                                    trace(outer_prod_sum(coul_adens * cb_bo, cb_bo)) +
                                    trace(outer_prod_sum(m_coul_bdens * ca_bo, ca_bo)) +
                                    trace(outer_prod_sum(m_coul_bdens * cb_bo, cb_bo));

    // Exchange energy:  Let K[P] be the Exchange matrix build from the density P
    // Again K[P] and P are block-diagonal
    //     energy =   tr( P^\alpha K[P^\alpha] ) + tr( P^\beta K[P^\beta] )
    const scalar_type energy_exchge =
          trace(outer_prod_sum(exchge_adens * ca_bo, ca_bo)) +
          trace(outer_prod_sum(m_exchge_bdens * cb_bo, cb_bo));

    // Scale energy appropriately and set it in map:
    // 0.5 because those are 2e term and we need to avoid double counting.
    base_type::m_energies[coul_adens.id()]   = 0.5 * coeff_coul * energy_coul;
    base_type::m_energies[exchge_adens.id()] = 0.5 * coeff_exchge * energy_exchge;
  }
}

template <typename StoredMatrix, typename BlockType>
void UnrestrictedIntegralOperator<StoredMatrix, BlockType>::assert_coefficient_structure(
      const coefficients_ptr_type& coeff_bf_ptr) const {
#ifdef DEBUG
  const size_t n_orbs = coeff_bf_ptr->n_vectors();
  const size_t n_bas  = m_coul_bdens.n_rows();
  assert_internal(n_orbs % 2 == 0);
  const size_t n_orbs_alpha = n_orbs / 2;

  auto is_zero = [](scalar_type v) { return v == 0.; };

  // alpha-beta block:
  for (size_t i = 0; i < n_orbs_alpha; ++i) {
    auto& vec = (*coeff_bf_ptr)[i];
    const bool all_zero =
          std::all_of(vec.begin() + n_bas, vec.begin() + 2 * n_bas, is_zero);
    assert_internal(all_zero);
  }

  // beta-alpha block:
  for (size_t i = n_orbs_alpha; i < n_orbs; ++i) {
    auto& vec = (*coeff_bf_ptr)[i];

    const bool all_zero = std::all_of(vec.begin(), vec.begin() + n_bas, is_zero);
    assert_internal(all_zero);
  }
#endif  // DEBUG
}

}  // namespace molsturm
