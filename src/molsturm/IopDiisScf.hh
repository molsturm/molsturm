#pragma once
#include "IntegralOperatorBase.hh"
#include "IopDiisScfKeys.hh"
#include "ScfErrorLibrary.hh"
#include <gscf/PulayDiisScf.hh>
#include <iostream>
#include <iterator>

// TODO A lot of duplicated code with the plain solver
// Maybe one should find a generalisation and merge them?
// Maybe some kind of output wrapper or convergence wrapper?

namespace molsturm {

template <typename ProblemMatrix, typename OverlapMatrix>
struct IopDiisScfState : public gscf::PulayDiisScfState<ProblemMatrix, OverlapMatrix> {
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

  // Norm of the most recent Pulay error
  real_type last_error_norm;

  IopDiisScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<real_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<real_type>::invalid},
          last_error_norm{linalgwrap::Constants<real_type>::invalid} {}
};

template <typename IntegralOperator, typename OverlapMatrix>
class IopDiisScf
      : public gscf::PulayDiisScf<IopDiisScfState<IntegralOperator, OverlapMatrix>> {
public:
  typedef IntegralOperator operator_type;
  typedef gscf::PulayDiisScf<IopDiisScfState<IntegralOperator, OverlapMatrix>> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  /** Construct scf from parameter map */
  IopDiisScf(const krims::ParameterMap& map) { update_control_params(map); }

  /** \name Iteration control */
  ///@{
  //! Maximal total energy change between two cycles.
  real_type max_tot_energy_change = 1e-8;

  //! Maximal 1e energy change between two cycles.
  real_type max_1e_energy_change = 1e-5;

  /** Check convergence by checking the maximal deviation of
   *  the last and previous eval pointers */
  bool is_converged(const state_type& state) const override {
    // We cannot be converged on the first iteration
    if (state.n_iter_count() <= 1) return false;

    // Check the most recent SCF error
    if (state.last_error_norm > base_type::max_error_norm) return false;

    // Total energy change
    const probmat_type& fock_bb = *state.problem_matrix_ptr();
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
          map.at(IopDiisScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change =
          map.at(IopDiisScfKeys::max_1e_energy_change, max_1e_energy_change);
  }
  ///@}

protected:
  matrix_type calculate_error(const state_type& s) const override {
    typedef ScfErrorLibrary<probmat_type> errorlib;
    return errorlib::pulay_error(s.overlap_matrix(), *s.eigenvectors_ptr(),
                                 *s.problem_matrix_ptr());
  }

  void before_iteration_step(state_type& s) const override {
    // Store the last step data away for use in the convergence check:
    s.last_step_1e_energy = s.problem_matrix_ptr()->energy_1e_terms();
    s.last_step_tot_energy =
          s.problem_matrix_ptr()->energy_2e_terms() + s.last_step_1e_energy;

    if (s.n_iter_count() == 1) {
      std::cout << "" << std::setw(5) << std::right << "iter" << std::setw(14) << "e1e"
                << std::setw(14) << "e2e" << std::setw(14) << "etot" << std::setw(14)
                << "scf_error" << std::endl;
    }
  }

  void after_iteration_step(state_type& s) const override {
    const probmat_type& fock_bb = *s.problem_matrix_ptr();
    auto n_iter = s.n_iter_count();

    // Compute the SCF Pulay error
    s.last_error_norm = norm_frobenius(s.errors.back());

    // scf_iter        e1e         e2e       etot        scf_error
    std::cout << " " << std::setw(4) << std::right << n_iter << std::setw(14)
              << fock_bb.energy_1e_terms() << std::setw(14) << fock_bb.energy_2e_terms()
              << std::setw(14) << fock_bb.energy_total() << std::setw(14)
              << s.last_error_norm << std::endl;
  }

  void on_converged(state_type& s) const override {
    auto& problem_matrix = *s.problem_matrix_ptr();
    auto n_iter = s.n_iter_count();

    // Indention unit:
    const std::string ind = "   ";

    std::cout << "Converged after " << n_iter << " iterations." << std::endl
              << std::endl
              << ind << "Hartree-Fock energies:" << std::endl;

    size_t longestfirst = 0;
    for (auto kv : problem_matrix.energies()) {
      longestfirst = std::max(kv.first.size(), longestfirst);
    }

    // Store original precision:
    linalgwrap::io::OstreamState outstate(std::cout);

    for (auto kv : problem_matrix.energies()) {
      std::cout << ind << ind << std::left << std::setw(longestfirst) << kv.first << " = "
                << kv.second << std::endl;
    }

    std::cout << ind << ind << std::left << std::setw(longestfirst) << "E_1e"
              << " = " << std::setprecision(10) << problem_matrix.energy_1e_terms()
              << std::endl
              << ind << ind << std::left << std::setw(longestfirst) << "E_2e"
              << " = " << std::setprecision(10) << problem_matrix.energy_2e_terms()
              << std::endl
              << std::endl
              << std::endl
              << ind << ind << std::left << std::setw(longestfirst) << "E_total"
              << " = " << std::setprecision(15) << problem_matrix.energy_total()
              << std::endl;
  }
};

}  // namespace molsturm
