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
class IopScfStateWrapper final : public InnerState {
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

  // Total energy of the previous step
  real_type last_step_tot_energy;

  // One electron energy of the previous step
  real_type last_step_1e_energy;

  IopScfStateWrapper(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<real_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<real_type>::invalid} {}
};

/** Wrapper around other SCF solvers, which IopScf uses to
 *  produce output and take care of a couple of things.
 */
template <typename InnerScf>
class IopScfWrapper final : public InnerScf {
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

  /** Check convergence by checking the maximal deviation of
   *  the last and previous eval pointers */
  bool is_converged(const state_type& state) const override;

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);

    max_tot_energy_change =
          map.at(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change = map.at(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    verbosity = map.at(IopScfKeys::verbosity, verbosity);
  }

  /** Get the current settings of all internal control parameters and
   *  update the ParameterMap accordingly.
   */
  void get_control_params(krims::ParameterMap& map) const {
    base_type::get_control_params(map);
    map.update(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    map.update(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    map.update(IopScfKeys::verbosity, verbosity);
  }
  ///@}

  IopScfWrapper(const krims::ParameterMap& map) : InnerScf(map) {
    update_control_params(map);
  }

 protected:
  matrix_type calculate_error(const state_type& s) const override;
  void before_iteration_step(state_type& s) const override;
  void after_iteration_step(state_type& s) const override;
  void on_converged(state_type& s) const override;
};

//
// -------------------------------------------------------
//
template <typename InnerScf>
bool IopScfWrapper<InnerScf>::is_converged(const state_type& state) const {
  // We cannot be converged on the first iteration
  if (state.n_iter() <= 1) return false;

  // Check the most recent SCF error
  if (state.last_error_norm > base_type::max_error_norm) return false;

  // Total energy change
  const probmat_type& fock_bb = state.problem_matrix();
  real_type tot_energy = fock_bb.energy_1e_terms() + fock_bb.energy_2e_terms();
  real_type tot_energy_diff = std::abs(state.last_step_tot_energy - tot_energy);
  if (tot_energy_diff > max_tot_energy_change) return false;

  // 1e energy change
  real_type e1_energy_diff =
        std::abs(state.last_step_1e_energy - fock_bb.energy_1e_terms());
  if (e1_energy_diff > max_1e_energy_change) return false;

  // TODO check convergence in density

  return true;
}

template <typename InnerScf>
typename IopScfWrapper<InnerScf>::matrix_type IopScfWrapper<InnerScf>::calculate_error(
      const state_type& s) const {
  // Forward pulay error as the error matrix for the DIIS and
  // other SCFs
  return ScfErrorLibrary<probmat_type>::pulay_error(
        s.overlap_matrix(), s.eigensolution().evectors(), s.problem_matrix());
}

template <typename InnerScf>
void IopScfWrapper<InnerScf>::before_iteration_step(state_type& s) const {
  // Store the last step data away for use in the convergence check:
  s.last_step_1e_energy = s.problem_matrix().energy_1e_terms();
  s.last_step_tot_energy = s.problem_matrix().energy_2e_terms() + s.last_step_1e_energy;

  if (s.n_iter() == 1) {
    if (have_common_bit(verbosity, ScfMsgType::IterationProcess)) {
      std::cout << std::setw(5) << std::right << "iter" << std::setw(12) << "e1e"
                << std::setw(12) << "e2e" << std::setw(12) << "etot" << std::setw(14)
                << "scf_error" << std::right << std::setw(12) << "n_eprob_it"
                << std::endl;
    }
  }
}

template <typename InnerScf>
void IopScfWrapper<InnerScf>::after_iteration_step(state_type& s) const {
  const probmat_type& fock_bb = s.problem_matrix();

  if (have_common_bit(verbosity, ScfMsgType::IterationProcess)) {
    // scf_iter        e1e         e2e       etot        scf_error
    std::cout << " " << std::setw(4) << std::right << s.n_iter() << std::setw(12)
              << fock_bb.energy_1e_terms() << std::setw(12) << fock_bb.energy_2e_terms()
              << std::setw(12) << fock_bb.energy_total() << std::setw(14)
              << s.last_error_norm << std::right << std::setw(12)
              << s.eigenproblem_stats().n_iter() << std::endl;
  }
}

template <typename InnerScf>
void IopScfWrapper<InnerScf>::on_converged(state_type& s) const {
  auto& fock_bb = s.problem_matrix();

  if (have_common_bit(verbosity, ScfMsgType::FinalSummary)) {
    // Indention unit:
    const std::string ind = "      ";

    std::cout << std::endl
              << "Converged after  " << std::endl
              << ind << "SCF iterations:    " << std::right << std::setw(7) << s.n_iter()
              << std::endl
              << ind << "operator applies:  " << std::right << std::setw(7)
              << s.n_mtx_applies() << std::endl
              << "with energies" << std::endl;

    size_t longestfirst = 0;
    for (auto kv : fock_bb.energies()) {
      longestfirst = std::max(kv.first.size(), longestfirst);
    }

    // Store original precision:
    linalgwrap::io::OstreamState outstate(std::cout);

    for (auto kv : fock_bb.energies()) {
      std::cout << ind << std::left << std::setw(longestfirst) << kv.first << " = "
                << kv.second << std::endl;
    }

    std::cout << ind << std::left << std::setw(longestfirst) << "E_1e"
              << " = " << std::setprecision(10) << fock_bb.energy_1e_terms() << std::endl
              << ind << std::left << std::setw(longestfirst) << "E_2e"
              << " = " << std::setprecision(10) << fock_bb.energy_2e_terms() << std::endl
              << std::endl
              << ind << std::left << std::setw(longestfirst) << "E_total"
              << " = " << std::setprecision(15) << fock_bb.energy_total() << std::endl;
  }
}

}  // namespace detail
}  // namespace molsturm
