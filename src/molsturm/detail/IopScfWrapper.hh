#pragma once
#include "molsturm/IntegralOperatorBase.hh"
#include "molsturm/IopScfKeys.hh"
#include "molsturm/ScfErrorLibrary.hh"
#include "molsturm/ScfMsgType.hh"
#include <iostream>
#include <iterator>

namespace molsturm {
namespace detail {

/** Wrapper around other SCF states, which adds things
 *  which are specific to IopScf
 */
template <typename InnerState>
class IopScfStateWrapper /*final*/ : public InnerState {
 public:
  // Use final here, because overwriting from the methods
  // makes no sense (we are the top level of the tree
  // implementing all the handlers)

  typedef InnerState base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::real_type real_type;

  static_assert(IsIntegralOperator<probmat_type>::value,
                "IopDiisScf only works sensibly with a proper IntegralOperator");

  // Total energy of the previous SCF iteration
  real_type last_step_tot_energy;

  // One electron energy of the previous SCF iteration
  real_type last_step_1e_energy;

  // Most recently encountered total energy change
  // (i.e. between the last iteration and the one before that)
  real_type last_tot_energy_change;

  // Most recently encountered one electron  energy change
  // (i.e. between the last iteration and the one before that)
  real_type last_1e_energy_change;

  IopScfStateWrapper(probmat_type probmat, const overlap_type& overlap_mat,
                     size_t n_iter_offset = 0)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<real_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<real_type>::invalid},
          m_n_iter_offset(n_iter_offset) {}

  /** Move the most recently obtained error values from another state */
  template <typename OtherState>
  void obtain_last_errors_from(const OtherState& s) {
    last_step_tot_energy = s.last_step_tot_energy;
    last_step_1e_energy = s.last_step_1e_energy;
    base_type::last_error_norm = s.last_error_norm;
    last_tot_energy_change = s.last_tot_energy_change;
    last_1e_energy_change = s.last_1e_energy_change;
  }

  size_t n_iter() const override { return m_n_iter_offset + base_type::n_iter(); }

 private:
  size_t m_n_iter_offset;
};

/** Wrapper around other SCF solvers, which IopScf uses to
 *  produce output and take care of a couple of things.
 */
template <typename InnerScf>
class IopScfWrapper /*final*/ : public InnerScf {
 public:
  // Use final here, because overwriting from the methods
  // makes no sense (we are the top level of the tree
  // implementing all the handlers)

  typedef InnerScf base_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::matrix_type matrix_type;
  /*
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::vector_type vector_type;
  */

  static_assert(
        std::is_same<state_type,
                     IopScfStateWrapper<typename state_type::base_type>>::value,
        "The state inside InnerScf needs to be an IopScfStateWrapper - wrapped state.");

  /** \name Iteration control */
  ///@{
  //! Maximal total energy change between two cycles.
  real_type max_tot_energy_change = 1e-8;

  //! Maximal 1e energy change between two cycles.
  real_type max_1e_energy_change = 1e-5;

  //! How verbose should the solver be.
  ScfMsgType verbosity = ScfMsgType::Silent;

  /** Run until this iteration count has been reached but no further.
   *
   * This is needed to switch to a different solver sensibly once
   * we have reached a certain iteration number. The special value
   * 0 indicates to never switch to another solver and finish
   * the SCF with the one which is currently running.
   **/
  size_t run_until_iter = 0;

  /** Check convergence */
  bool is_converged(const state_type& state) const override;

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);

    max_tot_energy_change =
          map.at(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change = map.at(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    verbosity = map.at(IopScfKeys::verbosity, verbosity);
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    map.update(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    map.update(IopScfKeys::verbosity, verbosity);
  }
  ///@}

  IopScfWrapper(const krims::GenMap& map) : InnerScf() { update_control_params(map); }

 protected:
  matrix_type calculate_error(const state_type& s) const override;
  void after_iteration_step(state_type& s) const override;
};

//
// -------------------------------------------------------
//
template <typename InnerScf>
bool IopScfWrapper<InnerScf>::is_converged(const state_type& state) const {
  // We cannot be converged on the first iteration
  if (state.n_iter() <= 1) return false;

  // If the run_until_iter iteration count has been reached, we consider
  // this solver done, i.e. converged.
  if (state.n_iter() >= run_until_iter) return true;

  if (state.last_error_norm > base_type::max_error_norm) return false;
  if (state.last_tot_energy_change > max_tot_energy_change) return false;
  if (state.last_1e_energy_change > max_1e_energy_change) return false;
  // TODO check convergence in density / coefficients
  // TODO code duplication with IopScf

  return true;
}

template <typename InnerScf>
typename IopScfWrapper<InnerScf>::matrix_type IopScfWrapper<InnerScf>::calculate_error(
      const state_type& s) const {
  // TODO: Note that this actually is an apply and should increase the apply count

  // Forward pulay error as the error matrix for the DIIS and other SCFs
  return ScfErrorLibrary<probmat_type>::pulay_error(s.problem_matrix(),
                                                    s.overlap_matrix());
}

template <typename InnerScf>
void IopScfWrapper<InnerScf>::after_iteration_step(state_type& s) const {
  const probmat_type& fock_bb = s.problem_matrix();

  // Update last step energy values and differences:
  s.last_1e_energy_change = std::abs(s.last_step_1e_energy - fock_bb.energy_1e_terms());
  s.last_tot_energy_change = std::abs(s.last_step_tot_energy - fock_bb.energy_total());
  s.last_step_1e_energy = s.problem_matrix().energy_1e_terms();
  s.last_step_tot_energy = s.problem_matrix().energy_total();

  if (have_common_bit(verbosity, ScfMsgType::IterationProcess)) {
    // scf_iter        e1e         e2e       etot        scf_error
    std::cout << " " << std::setw(4) << std::right << s.n_iter() << std::setw(12)
              << fock_bb.energy_1e_terms() << std::setw(12) << fock_bb.energy_2e_terms()
              << std::setw(12) << fock_bb.energy_total() << std::setw(14)
              << s.last_error_norm << std::right << std::setw(12)
              << s.eigenproblem_stats().n_iter() << std::endl;
  }
}

}  // namespace detail
}  // namespace molsturm
