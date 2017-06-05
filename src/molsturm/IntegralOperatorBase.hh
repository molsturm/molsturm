#pragma once
#include "IntegralTermContainer.hh"
#include "MolecularSystem.hh"
#include "energy_nuclear_repulsion.hh"
#include <gscf/FocklikeMatrix_i.hh>
#include <linalgwrap/LazyMatrix_i.hh>
#include <map>

namespace molsturm {

template <typename StoredMatrix>
class IntegralOperatorBase
      : public gscf::FocklikeMatrix_i<typename StoredMatrix::scalar_type> {
 public:
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;
  typedef typename stored_matrix_type::vector_type vector_type;

  /** Type of the coefficients and the pointer to the coefficients
   *  used inside the operator */
  typedef const linalgwrap::MultiVector<vector_type> coefficients_type;
  typedef std::shared_ptr<coefficients_type> coefficients_ptr_type;

  //! Type of a integral term
  typedef typename IntegralTermContainer<stored_matrix_type>::int_term_type int_term_type;

  /** \name Construct an integral operator */
  IntegralOperatorBase(IntegralTermContainer<StoredMatrix> integral_terms,
                       const MolecularSystem& system);

  // TODO implement general integral operator stuff in here

  /** Obtain the coefficients which were used to update the state of this operator,
   *  i.e. this returns the coefficients which have been most recently passed
   *  to this operator via update()
   */
  const coefficients_type& coefficients() const { return *m_coefficients_ptr; }

  /** Return the update key for the solver */
  const std::string& scf_update_key() const override final { return m_update_key; }

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   *
   * It expects the new coefficients under the parameter key
   * returned by scf_update_key()
   */
  void update(const krims::GenMap& map);

  /** Update the inner state:
   * Build the Fock matrix with the new coefficients
   */
  void update(const coefficients_ptr_type& coefficients);

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
  scalar_type energy_1e_terms() const override final;

  /** Return the sum of the energies of all two-electron terms */
  scalar_type energy_2e_terms() const override final {
    return m_energies.at(m_coul_adens.id()) + m_energies.at(m_exchge_adens.id());
  }

  /** Return the nuclear repulsion component */
  scalar_type energy_nuclear_repulsion() const { return m_nuclear_repulsion_energy; }

  /** Return the total energy
   *
   * This means the *electronic* and the *nuclear* component. */
  scalar_type energy_total() const {
    return energy_1e_terms() + energy_2e_terms() + energy_nuclear_repulsion();
  }

  /** Return the 1e terms
   *
   * \note that only the alpha-alpha or the beta-beta block is returned, i.e.
   *       the dimensionality of the returned matrix object is identical
   *       to the dimensionality of terms_alpha and terms_beta.
   * */
  virtual std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
  terms_1e() const;

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of alpha
   * spin. */
  virtual std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
  terms_alpha() const = 0;

  // TODO it would be really great to have a better interface
  // for the functions above and the one below.

  /** Return a map from the id strings of the integral terms to const
   * references to the lazy matrix objects, which represent the terms of beta
   * spin. */
  virtual std::map<gint::IntegralIdentifier, linalgwrap::LazyMatrixProduct<StoredMatrix>>
  terms_beta() const = 0;

 protected:
  /** Update the state of the lazy matrix operator terms and rebuild
   *  the operator m_operator from them.    */
  virtual void update_operator(const coefficients_ptr_type& coeff_bf_ptr) = 0;

  /** Recompute the energies from the current state of the terms
   * (First call update_operator on an update! )    */
  virtual void update_energies(const coefficients_ptr_type& coeff_bf_ptr) = 0;

  //

  /** Set of coefficients for the one electron integral terms */
  std::vector<scalar_type> m_coeff_1e;

  /** The one electron terms */
  std::vector<int_term_type> m_terms_1e;

  /** The coefficient for the coulomb term */
  scalar_type m_coeff_coul;

  /** The alpha coulomb term
   * The coulombic part resulting from the alpha electron density
   * \note This is not the part which sits in the alpha-alpha block !
   */
  int_term_type m_coul_adens;

  /** The coefficient for the HF exchange term */
  scalar_type m_coeff_exchge;

  /** The alpha exchange term
   * The coulombic 2e part resulting from the alpha electron density
   * \note This is not the part which sits in the alpha-alpha block !
   */
  int_term_type m_exchge_adens;

  //

  /** The current electronic energies */
  std::map<gint::IntegralIdentifier, scalar_type> m_energies;

  /** The nuclear repulsion energy */
  scalar_type m_nuclear_repulsion_energy;

  //

  /** Number of alpha electrons */
  size_t m_n_alpha;

  /** Number of alpha electrons */
  size_t m_n_beta;

  /** The current coefficients the fock operator is initialised to */
  coefficients_ptr_type m_coefficients_ptr;

  //! Key used for updating the state.
  const std::string m_update_key = "evec_coefficients";
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is an IntegralOperator
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef scalar_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsIntegralOperator : public std::false_type {};

template <typename T>
struct IsIntegralOperator<T, krims::VoidType<typename T::stored_matrix_type>>
      : public std::is_base_of<IntegralOperatorBase<typename T::stored_matrix_type>, T> {
};
//@}

//
// -------------------------------------------------------------------
//

template <typename StoredMatrix>
IntegralOperatorBase<StoredMatrix>::IntegralOperatorBase(
      IntegralTermContainer<StoredMatrix> integral_terms, const MolecularSystem& system)
      : m_coeff_1e{std::move(integral_terms.coefficients_1e)},
        m_terms_1e{integral_terms.integral_terms_1e},
        m_coeff_coul{integral_terms.coefficient_coulomb},
        m_coul_adens{integral_terms.coulomb_term},
        m_coeff_exchge{integral_terms.coefficient_exchange},
        m_exchge_adens{integral_terms.exchange_term},
        m_energies{},
        m_nuclear_repulsion_energy{molsturm::energy_nuclear_repulsion(system.structure)},
        m_n_alpha{system.n_alpha},
        m_n_beta{system.n_beta} {
  using namespace linalgwrap;

  // Check that alpha is equal to beta
  assert_equal(system.n_beta, system.n_alpha);

  // Check that number of terms and number of coefficients agrees:
  assert_size(m_terms_1e.size(), m_coeff_1e.size());

#ifdef DEBUG
  // Check operator size and zero energy terms
  const size_t op_size = m_coul_adens.n_rows();
#endif

  auto itterm = std::begin(m_terms_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm) {
    assert_size(itterm->n_cols(), op_size);
    assert_size(itterm->n_rows(), op_size);

    m_energies.insert(std::make_pair(itterm->id(), Constants<scalar_type>::zero));
  }

  assert_size(m_coul_adens.n_cols(), op_size);
  assert_size(m_coul_adens.n_rows(), op_size);
  m_energies.insert(std::make_pair(m_coul_adens.id(), Constants<scalar_type>::zero));

  assert_size(m_exchge_adens.n_cols(), op_size);
  assert_size(m_exchge_adens.n_rows(), op_size);
  m_energies.insert(std::make_pair(m_exchge_adens.id(), Constants<scalar_type>::zero));

  m_nuclear_repulsion_energy = molsturm::energy_nuclear_repulsion(system.structure);
}

template <typename StoredMatrix>
void IntegralOperatorBase<StoredMatrix>::update(const krims::GenMap& map) {
  if (map.exists(m_update_key)) {
    auto coeff_bf_ptr =
          static_cast<coefficients_ptr_type>(map.at_ptr<coefficients_type>(m_update_key));
    update(coeff_bf_ptr);
  }
}

template <typename StoredMatrix>
void IntegralOperatorBase<StoredMatrix>::update(
      const coefficients_ptr_type& coeff_bf_ptr) {
  // Only update if we got new coefficients:
  if (coeff_bf_ptr == m_coefficients_ptr) return;
  m_coefficients_ptr = coeff_bf_ptr;

  update_operator(coeff_bf_ptr);
  update_energies(coeff_bf_ptr);
}

template <typename StoredMatrix>
typename IntegralOperatorBase<StoredMatrix>::scalar_type
IntegralOperatorBase<StoredMatrix>::energy_1e_terms() const {
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
IntegralOperatorBase<StoredMatrix>::terms_1e() const {
  using linalgwrap::LazyMatrixProduct;
  std::map<gint::IntegralIdentifier, LazyMatrixProduct<StoredMatrix>> ret;

  auto itterm = std::begin(m_terms_1e);
  auto itcoeff = std::begin(m_coeff_1e);
  for (; itterm != std::end(m_terms_1e); ++itterm, ++itcoeff) {
    ret.insert(std::make_pair(itterm->id(), (*itcoeff) * (*itterm)));
  }

  return ret;
}

}  // namespace molsturm
