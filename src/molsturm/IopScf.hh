#pragma once
#include "IntegralOperatorBase.hh"
#include "IopScfKeys.hh"
#include "ScfErrorLibrary.hh"
#include "ScfMsgType.hh"
#include <gscf/PulayDiisScf.hh>
#include <iostream>
#include <iterator>

namespace molsturm {

template <typename ProblemMatrix, typename OverlapMatrix>
struct IopScfState : public gscf::PulayDiisScfState<ProblemMatrix, OverlapMatrix> {
  typedef gscf::PulayDiisScfState<ProblemMatrix, OverlapMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  static_assert(IsIntegralOperator<probmat_type>::value,
                "IopDiisScf only works sensibly with a proper IntegralOperator");

  // Total energy of the previous step
  real_type last_step_tot_energy;

  // One electron energy of the previous step
  real_type last_step_1e_energy;

  /** The orbital energies of the SCF */
  const std::vector<scalar_type>& orbital_energies() const {
    return base_type::eigensolution().evalues();
  }

  /** The orbital coefficients of the SCF */
  const linalgwrap::MultiVector<vector_type>& orbital_coeff() const {
    return base_type::eigensolution().evectors();
  }

  IopScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<real_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<real_type>::invalid} {}
};
/** Scf Solver which should be used for integral operators.
 *
 * TODO
 * Right now this is exactly the PulayDiisScf. Later we probably want
 * something more clever, which first does some DIIS, but later
 * switches it off (or for something else) when all of a sudden
 * the errors get to small (and the linear system to be solved in the
 * DIIS becomes ill-conditioned
 */
template <typename IntegralOperator, typename OverlapMatrix>
class IopScf : public gscf::PulayDiisScf<IopScfState<IntegralOperator, OverlapMatrix>> {
public:
  typedef IntegralOperator operator_type;
  typedef gscf::PulayDiisScf<IopScfState<IntegralOperator, OverlapMatrix>> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** Construct scf from parameter map */
  IopScf(const krims::ParameterMap& map) { update_control_params(map); }

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
  bool is_converged(const state_type& state) const override {
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

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);

    max_tot_energy_change =
          map.at(IopScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change = map.at(IopScfKeys::max_1e_energy_change, max_1e_energy_change);
    verbosity = map.at(IopScfKeys::verbosity, verbosity);
  }
  ///@}

protected:
  matrix_type calculate_error(const state_type& s) const override {
    typedef ScfErrorLibrary<probmat_type> errorlib;
    return errorlib::pulay_error(s.overlap_matrix(), s.eigensolution().evectors(),
                                 s.problem_matrix());
  }

  void before_iteration_step(state_type& s) const override {
    // Store the last step data away for use in the convergence check:
    s.last_step_1e_energy = s.problem_matrix().energy_1e_terms();
    s.last_step_tot_energy = s.problem_matrix().energy_2e_terms() + s.last_step_1e_energy;

    if (s.n_iter() == 1) {
      if (common_bit(verbosity, ScfMsgType::IterationProcess)) {
        std::cout << std::setw(5) << std::right << "iter" << std::setw(12) << "e1e"
                  << std::setw(12) << "e2e" << std::setw(12) << "etot" << std::setw(14)
                  << "scf_error" << std::right << std::setw(12) << "n_eprob_it"
                  << std::endl;
      }
    }
  }

  void after_iteration_step(state_type& s) const override {
    const probmat_type& fock_bb = s.problem_matrix();

    if (common_bit(verbosity, ScfMsgType::IterationProcess)) {
      // scf_iter        e1e         e2e       etot        scf_error
      std::cout << " " << std::setw(4) << std::right << s.n_iter() << std::setw(12)
                << fock_bb.energy_1e_terms() << std::setw(12) << fock_bb.energy_2e_terms()
                << std::setw(12) << fock_bb.energy_total() << std::setw(14)
                << s.last_error_norm << std::right << std::setw(12)
                << s.eigenproblem_stats().n_iter() << std::endl;
    }
  }

  void on_converged(state_type& s) const override {
    auto& fock_bb = s.problem_matrix();

    if (common_bit(verbosity, ScfMsgType::FinalSummary)) {
      // Indention unit:
      const std::string ind = "      ";

      std::cout << std::endl
                << "Converged after  " << std::endl
                << ind << "SCF iterations:    " << std::right << std::setw(7)
                << s.n_iter() << std::endl
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
                << " = " << std::setprecision(10) << fock_bb.energy_1e_terms()
                << std::endl
                << ind << std::left << std::setw(longestfirst) << "E_2e"
                << " = " << std::setprecision(10) << fock_bb.energy_2e_terms()
                << std::endl
                << std::endl
                << ind << std::left << std::setw(longestfirst) << "E_total"
                << " = " << std::setprecision(15) << fock_bb.energy_total() << std::endl;
    }
  }
};

/** Try to find the self consistent field configuration for an integral operator
 *  and an overlap matrix
 */
template <typename IntegralOperator, typename OverlapMatrix>
typename IopScf<IntegralOperator, OverlapMatrix>::state_type run_scf(
      IntegralOperator iop, const OverlapMatrix& s, krims::ParameterMap map) {
  IopScf<IntegralOperator, OverlapMatrix> solver(map);
  return solver.solve(std::move(iop), s);
}

}  // namespace molsturm
