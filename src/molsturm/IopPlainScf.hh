#pragma once
#include "IntegralOperatorBase.hh"
#include "IopPlainScfKeys.hh"
#include "ScfErrorLibrary.hh"
#include <gscf/PlainScf.hh>
#include <iostream>
#include <iterator>
#include <linalgwrap/io.hh>

// TODO A lot of duplicated code with the plain solver
// Maybe one should find a generalisation and merge them?
// Maybe some kind of output wrapper or convergence wrapper?

namespace molsturm {

template <typename ProblemMatrix, typename OverlapMatrix>
struct IopPlainScfState : public gscf::PlainScfState<ProblemMatrix, OverlapMatrix> {
  typedef gscf::PlainScfState<ProblemMatrix, OverlapMatrix> base_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  static_assert(IsIntegralOperator<probmat_type>::value,
                "IopPlainScf only works sensibly with a proper IntegralOperator");

  // Total energy of the previous step
  real_type last_step_tot_energy;

  // One electron energy of the previous step
  real_type last_step_1e_energy;

  // Orbital eigenvalues of the previous step
  std::shared_ptr<std::vector<real_type>> last_step_eval_ptr;

  // Current Frobenius norm of the Pulay SCF error
  real_type last_error_norm;

  IopPlainScfState(probmat_type probmat, const overlap_type& overlap_mat)
        : base_type{std::move(probmat), overlap_mat},
          last_step_tot_energy{linalgwrap::Constants<real_type>::invalid},
          last_step_1e_energy{linalgwrap::Constants<real_type>::invalid},
          last_step_eval_ptr{nullptr},
          last_error_norm{linalgwrap::Constants<real_type>::invalid} {}
};

template <typename IntegralOperator, typename OverlapMatrix>
class IopPlainScf
      : public gscf::PlainScf<IopPlainScfState<IntegralOperator, OverlapMatrix>> {
public:
  typedef IntegralOperator operator_type;
  typedef gscf::PlainScf<IopPlainScfState<IntegralOperator, OverlapMatrix>> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::probmat_type probmat_type;
  typedef typename base_type::overlap_type overlap_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::vector_type vector_type;

  IopPlainScf(linalgwrap::io::DataWriter_i<scalar_type>& writer,
              const krims::ParameterMap& map)
        : m_writer(writer) {
    update_control_params(map);
  }

  IopPlainScf(linalgwrap::io::DataWriter_i<scalar_type>& writer) : m_writer(writer) {}

  /** \name Iteration control */
  ///@{
  /** \brief Maximal frobenius norm of the pulay error matrix for an scf step
   *
   * The error is estimated using the Pulay SCF error
   * in the atom centered basis */
  real_type max_error_norm = 5e-7;

  //! Maximal total energy change between two cycles.
  real_type max_tot_energy_change = 1e-8;

  //! Maximal 1e energy change between two cycles.
  real_type max_1e_energy_change = 1e-5;

  // Should orbital energies be printed?
  bool print_orbital_energies = false;

  /** Check convergence by checking the maximial deviation of
   *  the last and previous eval pointers */
  bool is_converged(const state_type& state) const override {
    // If no previous eval_ptr, then we are not converged:
    if (!state.last_step_eval_ptr) return false;

    // Check the SCF error
    if (state.last_error_norm > max_error_norm) return false;

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

    max_error_norm = map.at(IopPlainScfKeys::max_error_norm, max_error_norm);
    max_tot_energy_change =
          map.at(IopPlainScfKeys::max_tot_energy_change, max_tot_energy_change);
    max_1e_energy_change =
          map.at(IopPlainScfKeys::max_1e_energy_change, max_1e_energy_change);
    print_orbital_energies =
          map.at(IopPlainScfKeys::print_orbital_energies, print_orbital_energies);
  }
  ///@}

protected:
  void before_iteration_step(state_type& s) const override {
    // Store the last step data away for use in the convergence check:
    s.last_step_eval_ptr = s.eigenvalues_ptr();
    s.last_step_1e_energy = s.problem_matrix_ptr()->energy_1e_terms();
    s.last_step_tot_energy =
          s.problem_matrix_ptr()->energy_2e_terms() + s.last_step_1e_energy;

    if (s.n_iter_count() == 1) {
      std::cout << " " << std::setw(5) << std::right << "iter" << std::setw(14) << "e1e"
                << std::setw(14) << "e2e" << std::setw(14) << "etot" << std::setw(14)
                << "scf_error" << std::endl;
    }
  }

  void on_update_eigenpairs(state_type& s) const override {
    assert_dbg(m_writer, krims::ExcIO());
    m_writer.write("evals" + std::to_string(s.n_iter_count()),
                   linalgwrap::make_as_multivector<vector_type>(*s.eigenvalues_ptr()));
    m_writer.write("evecs" + std::to_string(s.n_iter_count()), *s.eigenvectors_ptr());

    // Print orbital evals:
    if (print_orbital_energies) {
      const std::string indent = "         ";
      std::cout << indent << "New orbital eigenvalues: " << std::endl;
      std::cout << indent;
      std::ostream_iterator<scalar_type> out_it(std::cout, " ");
      std::copy(s.eigenvalues_ptr()->begin(), s.eigenvalues_ptr()->end(), out_it);
      std::cout << std::endl;
    }

    assert_dbg(m_writer, krims::ExcIO());
  }

  void on_update_problem_matrix(state_type& s) const override {
    const probmat_type& fock_bb = *s.problem_matrix_ptr();
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
  }

  void after_iteration_step(state_type& s) const override {
    const probmat_type& fock_bb = *s.problem_matrix_ptr();
    auto n_iter = s.n_iter_count();

    // Compute the SCF Pulay error
    typedef ScfErrorLibrary<probmat_type> errorlib;
    s.last_error_norm = norm_frobenius(
          errorlib::pulay_error(s.overlap_matrix(), *s.eigenvectors_ptr(), fock_bb));

    // scf_iter        e1e         e2e       etot        scf_error
    std::cout << "" << std::setw(4) << std::right << n_iter << std::setw(14)
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

private:
  linalgwrap::io::DataWriter_i<scalar_type>& m_writer;
};

}  // namespace molsturm
