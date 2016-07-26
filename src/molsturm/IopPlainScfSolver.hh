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

  // Current Frobenius norm of the Pulay SCF error
  scalar_type pulay_error_norm;

  IopPlainScfSolverState(probmat_type probmat, const matrix_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<scalar_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<scalar_type>::invalid},
          last_step_eval_ptr{nullptr},
          pulay_error_norm{linalgwrap::Constants<scalar_type>::invalid} {}
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

  /** \brief Maximal frobenius norm of the pulay error matrix for an scf step
   *
   * The error is estimated using the Pulay SCF error
   * in the atom centered basis */
  scalar_type max_pulay_error_norm;

  // Maximal total energy change between two cycles.
  scalar_type max_tot_energy_change;

  // Maximal 1e energy change between two cycles.
  scalar_type max_1e_energy_change;

  IopPlainScfSolverControl()
        : base_type{},
          max_pulay_error_norm{5e-7},
          max_tot_energy_change{1e-8},
          max_1e_energy_change{1e-5} {}

  /** Check convergence by checking the maximial deviation of
   *  the last and previous eval pointers */
  bool is_converged(const scf_state_type& state) const override {
    // If no previous eval_ptr, then we are not converged:
    if (!state.last_step_eval_ptr) return false;

    // Check the SCF error
    if (state.pulay_error_norm > max_pulay_error_norm) return false;

    // Total energy change
    const probmat_type& fock_bb = *state.problem_matrix_ptr();
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
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

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
      std::cout << " " << std::setw(5) << std::right << "iter" << std::setw(14)
                << "e1e" << std::setw(14) << "e2e" << std::setw(14) << "etot"
                << std::setw(14) << "scf_error" << std::endl;
    }
  }

  void on_update_eigenpairs(scf_state_type& s) const override {
    std::cout << "New orbital eigenvalues: " << std::endl;

    assert_dbg(m_writer, linalgwrap::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   *s.eigenvalues_ptr());
    m_writer.write("evecs" + std::to_string(s.n_iter_count()),
                   *s.eigenvectors_ptr());

    // Print orbital evals:
    std::ostream_iterator<scalar_type> out_it(std::cout, " ");
    std::copy(s.eigenvalues_ptr()->begin(), s.eigenvalues_ptr()->end(), out_it);
    std::cout << std::endl;
  }

  void on_update_problem_matrix(scf_state_type& s) const override {
    const probmat_type& fock_bb = *s.problem_matrix_ptr();
    const matrix_type& coeff_bf = *s.eigenvectors_ptr();
    const matrix_type& overlap_bb = s.overlap_matrix();
    auto n_iter = s.n_iter_count();
    std::string itstr = std::to_string(n_iter);

    // Write current findings:
    for (auto kv : fock_bb.terms_alpha()) {
      // Normalise the label: The id of the term may contain funny symbols
      std::string lala = m_writer.normalise_label(kv.first + "a" + itstr);
      m_writer.write(lala, kv.second);
    }
    for (auto kv : fock_bb.terms_beta()) {
      // Normalise the label: The id of the term may contain funny symbols
      std::string lalb = m_writer.normalise_label(kv.first + "b" + itstr);
      m_writer.write(lalb, kv.second);
    }
    m_writer.write("fock" + itstr, fock_bb);

    // Compute the SCF Pulay error
    typedef ScfErrorLibrary<probmat_type> errorlib;
    s.pulay_error_norm =
          errorlib::pulay_error(overlap_bb, coeff_bf, fock_bb).norm_frobenius();

    // scf_iter        e1e         e2e       etot        scf_error
    std::cout << " " << std::setw(4) << std::right << n_iter << std::setw(14)
              << fock_bb.energy_1e_terms() << std::setw(14)
              << fock_bb.energy_2e_terms() << std::setw(14)
              << fock_bb.energy_total() << std::setw(14) << s.pulay_error_norm
              << std::endl;
  }

  void on_converged(scf_state_type& s) const override {
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
    // TODO use the tool i wrote for linalgwrap here!
    std::streamsize prec = std::cout.precision();

    for (auto kv : problem_matrix.energies()) {
      std::cout << ind << ind << std::left << std::setw(longestfirst)
                << kv.first << " = " << kv.second << std::endl;
    }

    std::cout << ind << ind << std::left << std::setw(longestfirst) << "E_1e"
              << " = " << std::setprecision(10)
              << problem_matrix.energy_1e_terms() << std::endl
              << ind << ind << std::left << std::setw(longestfirst) << "E_2e"
              << " = " << std::setprecision(10)
              << problem_matrix.energy_2e_terms() << std::endl
              << std::endl
              << std::endl
              << ind << ind << std::left << std::setw(longestfirst) << "E_total"
              << " = " << std::setprecision(14) << problem_matrix.energy_total()
              << std::endl;

    std::cout.precision(prec);
  }

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

}  // namespace molsturm
