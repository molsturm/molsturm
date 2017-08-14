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
#include "detail/IopScfWrapper.hh"
#include <gscf/PlainScf.hh>
#include <gscf/PulayDiisScf.hh>
#include <gscf/TruncatedOptDampScf.hh>

namespace molsturm {

// TODO Maybe a common base between this guy and lazyten's EigensystemSolver is
// sensible?

template <typename ProblemMatrix, typename OverlapMatrix>
struct IopScfState final : public gscf::ScfStateBase<ProblemMatrix, OverlapMatrix> {
  typedef gscf::ScfStateBase<ProblemMatrix, OverlapMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** The orbital energies of the SCF */
  const std::vector<scalar_type>& orbital_energies() const {
    return base_type::eigensolution().evalues();
  }

  /** The orbital coefficients of the SCF */
  const lazyten::MultiVector<vector_type>& orbital_coeff() const {
    return base_type::eigensolution().evectors();
  }

  IopScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat}, m_n_iter(0), m_n_mtx_applies(0) {}

  /** Transfer data from an arbitrary other state by copying
   * it into this one. This is needed, since we need to get the
   * data back here once the inner solvers are done with
   * (partially) solving the SCF.
   *
   * \note Do not use this to setup a guess state.
   * For this case the function obtain_guess_from exists.
   **/
  template <typename State>
  void push_intermediate_results(State&& other_state) {
    // Copy things from ScfStateBase
    // TODO This is all but good. The DiagMat template argument really has to
    //      disappear (or we need yet a layer ... a common base independent of
    //      DiagMat.
    base_type::last_error_norm    = other_state.last_error_norm;
    base_type::problem_matrix_ptr = other_state.problem_matrix_ptr;
    base_type::push_new_eigensolution(other_state.previous_eigensolution(),
                                      other_state.eigenproblem_stats());
    base_type::push_new_eigensolution(other_state.eigensolution(),
                                      other_state.eigenproblem_stats());

    last_tot_energy_change = other_state.last_tot_energy_change;
    last_1e_energy_change  = other_state.last_1e_energy_change;

    last_step_tot_energy = other_state.last_step_tot_energy;
    last_step_1e_energy  = other_state.last_step_1e_energy;

    // The following is correct, since we pass the current state number on to the
    // IopScfStateWrapper as well, such that it accumulates to give the correct
    // actual iteration count.
    m_n_iter = other_state.n_iter();
    m_n_mtx_applies += other_state.n_mtx_applies();

    static_cast<lazyten::SolverStateBase&>(*this) =
          static_cast<lazyten::SolverStateBase&&>(other_state);
  }

  size_t n_iter() const override { return m_n_iter; }
  size_t n_mtx_applies() const override { return m_n_mtx_applies; }

  // Total energy of the previous SCF iteration
  real_type last_step_tot_energy;

  // One electron energy of the previous SCF iteration
  real_type last_step_1e_energy;

  // Most recently encountered total energy change
  // (i.e. between the last inner scf iteration and the
  //  one before that)
  real_type last_tot_energy_change;

  // Most recently encountered one electron  energy change
  // (i.e. between the last inner scf iteration and the
  //  one before that)
  real_type last_1e_energy_change;

 private:
  size_t m_n_iter;
  size_t m_n_mtx_applies;
};

/** Scf Solver which should be used for integral operators.
 */
template <typename IntegralOperator, typename OverlapMatrix>
class IopScf final : public gscf::ScfBase<IopScfState<IntegralOperator, OverlapMatrix>> {
 public:
  typedef IntegralOperator operator_type;
  typedef gscf::ScfBase<IopScfState<IntegralOperator, OverlapMatrix>> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** Construct scf from parameter map */
  IopScf(const krims::GenMap& map) { update_control_params(map); }

  /** \name Iteration control */
  ///@{
  //! Maximal total energy change between two cycles.
  real_type max_tot_energy_change = 1e-8;

  //! Maximal 1e energy change between two cycles.
  real_type max_1e_energy_change = 1e-5;

  //! How verbose should the solver be.
  ScfMsgType verbosity = ScfMsgType::Silent;

  /** Print the progress during the solve
   *  (Type: bool)
   *
   *  Will add ScfMsgType::IterationProcess to the verbosity parameter.
   */
  static const std::string print_progress;

  // TODO introduce method parameter!

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);

    max_tot_energy_change =
          map.at(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change = map.at(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    verbosity            = map.at(IopScfKeys::verbosity, verbosity);

    // For transition ... later remove the verbosity flag.
    if (map.at(IopScfKeys::print_iterations, false)) {
      verbosity |= ScfMsgType::IterationProcess;
    }

    // Copy the map to the internal storage such that we
    // can pass it on to the actual eigensolvers.
    m_inner_params = map;
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    map.update(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    map.update(IopScfKeys::verbosity, verbosity);

    // For transition ... later remove the verbosity flag.
    map.update(IopScfKeys::print_iterations,
               have_common_bit(verbosity, ScfMsgType::IterationProcess));
  }
  ///@}

  /** Solve the state. Uses different SCF algorithms depending
   *  on certain properties of the iteration process */
  void solve_state(state_type& state) const override;

 protected:
  bool is_converged(const state_type& state) const override;
  void on_converged(state_type& s) const override;

 private:
  /** Types of the wrapped solvers */
  //@{
  typedef detail::IopScfWrapper<gscf::TruncatedOptDampScf<detail::IopScfStateWrapper<
        gscf::TruncatedOptDampScfState<IntegralOperator, OverlapMatrix>>>>
        TruncODASolver;
  typedef detail::IopScfWrapper<gscf::PulayDiisScf<detail::IopScfStateWrapper<
        gscf::PulayDiisScfState<IntegralOperator, OverlapMatrix>>>>
        DiisSolver;
  typedef detail::IopScfWrapper<gscf::PlainScf<
        detail::IopScfStateWrapper<gscf::PlainScfState<IntegralOperator, OverlapMatrix>>>>
        PlainSolver;
  //@}

  /** Solve until a provided maximal error norm is reached
   *
   *  Sets the error values as follows:
   *  - max_error_norm            the provided value
   *  - max_tot_energy_change     the provided value
   *  - max_1e_energy_change      100 times the provided value
   *  Of course: If the user wants only a less accurate solve, than
   *  this setting takes preference.
   *
   *  \param s           State type (guess on input, result on output)
   *  \param max_error   Maximum error to solve for. 0 means same
   *                     accuracy as user asked for via GenMap
   *                     or setting class arguments.
   *  \param run_until_iter    Maximum iteration count until which to run this solver.
   *                     I.e. if max_iter == 5 and the current
   *                     iteration to come is the 4th, then this solver
   *                     will be run at most for 2 steps (step 4 and 5).
   *                     The value of all is special and indicates, that
   *                     the solver will run until convergence in the error
   *                     has been reached.
   */
  template <typename WrappedSolver>
  void solve_up_to(state_type& s, real_type error,
                   size_t run_until_iter = lazyten::Constants<size_t>::all) const;

  /** Cache of the parameters which will be passed to the actual SCF.
   *
   * A mixture of the original parameters, which the user passed to us
   * and the current state which is reflected in the member variables
   * in this class and the subclasses. Should be updated with
   * get_control_params *before* the actual inner eigensolver invocation.
   */
  mutable krims::GenMap m_inner_params;
};

/** Try to find the self consistent field configuration for an integral operator
 *  and an overlap matrix
 */
template <typename IntegralOperator, typename OverlapMatrix>
typename IopScf<IntegralOperator, OverlapMatrix>::state_type run_scf(
      IntegralOperator iop, const OverlapMatrix& s, const krims::GenMap& map) {
  typedef IopScf<IntegralOperator, OverlapMatrix> scf;
  return scf{map}.solve(std::move(iop), s);
}

/** Try to find the self consistent field configuration for an integral operator
 *  and an overlap matrix given a guess solution.
 */
template <typename IntegralOperator, typename OverlapMatrix>
typename IopScf<IntegralOperator, OverlapMatrix>::state_type run_scf(
      IntegralOperator iop, const OverlapMatrix& s,
      typename IopScf<IntegralOperator, OverlapMatrix>::state_type::esoln_type
            guess_solution,
      const krims::GenMap& map);

//
// -----------------------------------------------------------
//

template <typename IntegralOperator, typename OverlapMatrix>
template <typename WrappedSolver>
void IopScf<IntegralOperator, OverlapMatrix>::solve_up_to(state_type& state,
                                                          real_type error,
                                                          size_t run_until_iter) const {
  // Setup inner solver (max => User settings of accuracy take preference)
  WrappedSolver inner_solver{m_inner_params};
  inner_solver.max_error_norm        = std::max(base_type::max_error_norm, error);
  inner_solver.max_1e_energy_change  = std::max(100. * error, max_1e_energy_change);
  inner_solver.max_tot_energy_change = std::max(error, max_tot_energy_change);
  inner_solver.run_until_iter        = run_until_iter;

  // Setup inner state
  typename WrappedSolver::state_type inner_state{state.problem_matrix(),
                                                 state.overlap_matrix(), state.n_iter()};

  if (state.n_iter() > 1) {
    inner_state.obtain_last_errors_from(state);

    // TODO Move state here -> gets rid of a copy!
    inner_state.obtain_guess_from(state);
  }

  // Do solve
  try {
    inner_solver.solve_state(inner_state);
    state.push_intermediate_results(std::move(inner_state));
  } catch (lazyten::SolverException& e) {
    // On exception still update the state reference
    state.push_intermediate_results(std::move(inner_state));
    throw;
  }
}

template <typename IntegralOperator, typename OverlapMatrix>
void IopScf<IntegralOperator, OverlapMatrix>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  // Make sure that control parameters are up to date:
  get_control_params(m_inner_params);

  const bool print_progress = have_common_bit(verbosity, ScfMsgType::IterationProcess);

  // Print table header if user wishes an iteration progress:
  if (print_progress) {
    std::cout.unsetf(std::ios::floatfield);  // Unset floatfield flags
    std::cout << std::setprecision(7) << std::setw(5) << std::right << "iter"
              << std::setw(12) << "e1e" << std::setw(12) << "e2e" << std::setw(12)
              << "etot" << std::setw(16) << "scf_error" << std::right << std::setw(12)
              << "n_eprob_it" << std::endl;
  }

  // TODO make this configurable
  const real_type diis_startup_error_norm = 0.25;
  const size_t diis_startup_iter          = 12;

  {  // truncated ODA
    solve_up_to<TruncODASolver>(state, diis_startup_error_norm, diis_startup_iter);
    if (base_type::convergence_reached(state)) return;
  }

  // TODO make this configurable
  //      Idea: have a mapping
  //         { {accuracy, method},
  //           {accuracy2, method2},
  //         }
  const real_type diis_limit_max_error_norm = 5e-7;

  {  // DIIS
    if (print_progress) {
      // TODO It would be nice to have a function to do this printing here
      std::cout << "                  ****    Turning on DIIS   ****" << std::endl;
    }
    solve_up_to<DiisSolver>(state, diis_limit_max_error_norm);
    if (base_type::convergence_reached(state)) return;
    if (print_progress) {
      std::cout << "                  ****  Switching off DIIS  ****" << std::endl;
    }
  }

  // TODO From this point we ideally want to do a SOSCF

  // First experiments seem to suggest that the tODA is not good
  // to get to extremely high accuracy due to numerical errors in the
  // computation of the traces. Hence we cap it at an error of 1e-7
  const real_type toda_limit_max_error_norm = 1e-7;

  {  // truncated ODA SCF
    solve_up_to<TruncODASolver>(state, toda_limit_max_error_norm);
    if (base_type::convergence_reached(state)) return;
  }

  {  // Plain SCF
    if (print_progress) {
      std::cout << "                  **** Removing any damping ****" << std::endl;
    }
    solve_up_to<PlainSolver>(state, 0.);
    if (base_type::convergence_reached(state)) return;
  }

  // Cannot happen since last solver should have converged or thrown
  assert_internal(false);
}

template <typename IntegralOperator, typename OverlapMatrix>
bool IopScf<IntegralOperator, OverlapMatrix>::is_converged(
      const state_type& state) const {
  if (state.last_error_norm > base_type::max_error_norm) return false;
  if (state.last_tot_energy_change > max_tot_energy_change) return false;
  if (state.last_1e_energy_change > max_1e_energy_change) return false;
  // TODO check convergence in density / coefficients
  // TODO code duplication with IopScfWrapper

  return true;
}

template <typename IntegralOperator, typename OverlapMatrix>
void IopScf<IntegralOperator, OverlapMatrix>::on_converged(state_type& s) const {
  auto& fock_bb = s.problem_matrix();

  if (have_common_bit(verbosity, ScfMsgType::FinalSummary)) {
    // Indention unit:
    const std::string ind = "      ";
    const std::string nuc_rep_label("nuclear repulsion");

    std::cout.unsetf(std::ios::floatfield);  // Unset floatfield flags
    std::cout << std::endl
              << "Converged after  " << std::endl
              << ind << "SCF iterations:    " << std::right << std::setw(7) << s.n_iter()
              << std::endl
              << ind << "operator applies:  " << std::right << std::setw(7)
              << s.n_mtx_applies() << std::endl
              << "with energies" << std::endl;

    // For the virial ratio we need accumulated kinetic and potential energies.
    real_type kinetic_energy   = 0;
    real_type potential_energy = fock_bb.energy_nuclear_repulsion();

    size_t friendly_label_size = nuc_rep_label.size();
    for (const auto& kv : fock_bb.energies()) {
      const auto& id = kv.first;

      friendly_label_size =
            std::max(id.integral_friendly_name().size(), friendly_label_size);

      if (id.integral_type() == gint::IntegralType::kinetic) {
        kinetic_energy += kv.second;
      } else {
        potential_energy += kv.second;
      }
    }

    // Store original precision:
    lazyten::io::OstreamState outstate(std::cout);

    // Print energy terms
    std::cout << ind << std::left << std::setw(friendly_label_size) << nuc_rep_label
              << " = " << std::setprecision(10) << fock_bb.energy_nuclear_repulsion()
              << '\n';

    for (const auto& kv : fock_bb.energies()) {
      std::cout << ind << std::left << std::setw(friendly_label_size)
                << kv.first.integral_friendly_name() << " = " << kv.second << '\n';
    }

    // Print 1e, 2e energies and total electronic energy
    const real_type E_electronic = fock_bb.energy_1e_terms() + fock_bb.energy_2e_terms();
    std::cout << '\n'
              << ind << std::left << std::setw(friendly_label_size) << "E_1e"
              << " = " << std::setprecision(10) << fock_bb.energy_1e_terms() << '\n'
              << ind << std::left << std::setw(friendly_label_size) << "E_2e"
              << " = " << std::setprecision(10) << fock_bb.energy_2e_terms() << '\n'
              << ind << std::left << std::setw(friendly_label_size) << "E electronic"
              << " = " << std::setprecision(10) << E_electronic << '\n';

    // Print kinetic and potential and virial ratio
    const real_type virial = -potential_energy / kinetic_energy;
    std::cout << '\n'
              << ind << std::left << std::setw(friendly_label_size) << "E_pot"
              << " = " << std::setprecision(10) << potential_energy << '\n'
              << ind << std::left << std::setw(friendly_label_size) << "E_kin"
              << " = " << std::setprecision(10) << kinetic_energy << '\n'
              << ind << std::left << std::setw(friendly_label_size) << "virial ratio"
              << " = " << std::setprecision(10) << virial << '\n'
              << '\n';

    std::cout << ind << std::left << std::setw(friendly_label_size) << "E_total"
              << " = " << std::setprecision(15) << fock_bb.energy_total() << std::endl;
  }
}

//
// -------------------------------------------------------------------
//

template <typename IntegralOperator, typename OverlapMatrix>
typename IopScf<IntegralOperator, OverlapMatrix>::state_type run_scf(
      IntegralOperator iop, const OverlapMatrix& s,
      typename IopScf<IntegralOperator, OverlapMatrix>::state_type::esoln_type
            guess_solution,
      const krims::GenMap& map) {
  typedef IopScf<IntegralOperator, OverlapMatrix> scf;

  IopScfState<IntegralOperator, OverlapMatrix> state(std::move(iop), s);
  state.obtain_guess_from(std::move(guess_solution));

  scf{map}.solve_state(state);
  return state;
}

}  // namespace molsturm
