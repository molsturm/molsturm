#pragma once
#include "ScfErrorLibrary.hh"
#include <gscf/PlainScf.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

namespace molsturm {

template <typename ScfTraits>
struct IopPlainScfSolverState : public gscf::PlainScfState<ScfTraits> {
  typedef gscf::PlainScfState<ScfTraits> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  static_assert(
        IsIntegralOperator<probmat_type>::value,
        "IopPlainScfSolver only works sensibly with a proper IntegralOperator");

  // Total energy of the previous step
  scalar_type last_step_tot_energy;

  // One electron energy of the previous step
  scalar_type last_step_1e_energy;

  // Previous orbital coefficients (evals of the problem matrix)
  std::shared_ptr<vector_type> last_step_eval_ptr;

  IopPlainScfSolverState(probmat_type probmat, const matrix_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<scalar_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<scalar_type>::invalid},
          last_step_eval_ptr{nullptr} {}
};

template <typename ScfState>
struct IopPlainScfSolverControl : public gscf::PlainScfControl<ScfState> {
  typedef gscf::PlainScfControl<ScfState> base_type;
  typedef typename base_type::scf_state_type scf_state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  static_assert(
        IsIntegralOperator<probmat_type>::value,
        "IopPlainScfSolver only works sensibly with a proper IntegralOperator");

  /** \brief Maximal l1 norm of the pulay error matrix for an scf step
   *
   *
   * The error is estimated using the Pulay SCF error
   * in the atom centered basis */
  scalar_type max_l1_pulay_error;

  // Maximal total energy change between two cycles.
  scalar_type max_tot_energy_change;

  // Maximal 1e energy change between two cycles.
  scalar_type max_1e_energy_change;

  IopPlainScfSolverControl()
        : base_type{},
          max_l1_pulay_error{5e-7},
          max_tot_energy_change{1e-8},
          max_1e_energy_change{1e-5} {}

  /** Check convergence by checking the maximial deviation of
   *  the last and previous eval pointers */
  bool is_converged(const scf_state_type& state) const override {
    // If no previous eval_ptr, then we are not converged:
    if (!state.last_step_eval_ptr) return false;

    // Extract references to relevant quantities:
    const probmat_type& fock_bb = *state.problem_matrix_ptr();
    const matrix_type& coeff_bf = *state.eigenvectors_ptr();
    const matrix_type& overlap_bb = state.overlap_matrix();

    // Check the SCF error
    typedef ScfErrorLibrary<probmat_type> errorlib;
    scalar_type l1_pulay_error =
          errorlib::pulay_error(overlap_bb, coeff_bf, fock_bb).norm_l1();
    if (l1_pulay_error > max_l1_pulay_error) return false;

    // Total energy change
    scalar_type tot_energy =
          fock_bb.energy_1e_terms() + fock_bb.energy_2e_terms();
    scalar_type tot_energy_diff =
          std::abs(state.last_step_tot_energy - tot_energy);
    if (tot_energy_diff > max_tot_energy_change) return false;

    // 1e energy change
    scalar_type e1_energy_diff =
          std::abs(state.last_step_1e_energy - fock_bb.energy_1e_terms());
    if (e1_energy_diff > max_1e_energy_change) return false;

    // TODO check convergence in density

    return true;
  }
};

template <typename IntegralOperator>
class IopPlainScfSolver
      : public gscf::PlainScf<IntegralOperator,
                              IopPlainScfSolverState<IntegralOperator>,
                              IopPlainScfSolverControl<
                                    IopPlainScfSolverState<IntegralOperator>>> {
public:
  typedef IntegralOperator operator_type;
  typedef gscf::PlainScf<
        IntegralOperator, IopPlainScfSolverState<IntegralOperator>,
        IopPlainScfSolverControl<IopPlainScfSolverState<IntegralOperator>>>
        base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::scf_state_type scf_state_type;

  IopPlainScfSolver(linalgwrap::io::DataWriter_i<scalar_type>& writer)
        : m_writer(writer) {}

protected:
  void before_iteration_step(scf_state_type& s) const override {
    // Store the last step data away for use in the convergence check:
    s.last_step_eval_ptr = s.eigenvalues_ptr();
    s.last_step_1e_energy = s.problem_matrix_ptr()->energy_1e_terms();
    s.last_step_tot_energy =
          s.problem_matrix_ptr()->energy_2e_terms() + s.last_step_1e_energy;

    if (s.n_iter_count() == 1) {
      std::cout
            << " scf_iter        e1e         e2e       etot        scf_error";
    }
  }

  void on_update_eigenpairs(scf_state_type& s) const override {
    std::cout << "       New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, linalgwrap::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   *s.eigenvalues_ptr());
    m_writer.write("evecs" + std::to_string(s.n_iter_count()),
                   *s.eigenvectors_ptr());

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::cout << "        ";
    std::copy(s.eigenvalues_ptr()->begin(), s.eigenvalues_ptr()->end(), out_it);
    std::cout << std::endl;
  }

  void on_update_problem_matrix(scf_state_type& s) const override {
    auto& problem_matrix = *s.problem_matrix_ptr();
    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    for (auto kv : problem_matrix.terms_alpha()) {
      m_writer.write(kv.first + "a" + itstr, kv.second);
    }
    for (auto kv : problem_matrix.terms_beta()) {
      m_writer.write(kv.first + "b" + itstr, kv.second);
    }
    m_writer.write("fock" + itstr, problem_matrix);

    // Calculate error
    scalar_type error = 0.0;

    // scf_iter        e1e         e2e       etot        scf_error
    std::cout << n_iter << "  " << problem_matrix.energy_1e_terms() << "  "
              << problem_matrix.energy_2e_terms() << "    "
              << problem_matrix.energy_total() << "    " << error << std::endl;
  }

  void on_converged(scf_state_type& s) const override {
    auto& problem_matrix = *s.problem_matrix_ptr();
    auto n_iter = s.n_iter_count();

    std::cout << std::endl
              << "Converged after " << n_iter << " iterations." << std::endl
              << std::endl
              << "   HF energies:" << std::endl;
    for (auto kv : problem_matrix.energies()) {
      std::cout << "     " << kv.first << " = " << kv.second << std::endl;
    }
    std::cout << "     E_1e     = " << problem_matrix.energy_1e_terms()
              << std::endl
              << "     E_2e     = " << problem_matrix.energy_2e_terms()
              << std::endl
              << std::endl
              << "     E_total  = " << problem_matrix.energy_total()
              << std::endl;
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

}  // namespace molsturm
