#include "parse_args.hh"

#include <gint/IntegralLookup.hh>
#include <gint/version.hh>
#include <gscf/version.hh>
#include <iostream>
#include <linalgwrap/io.hh>
#include <linalgwrap/version.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/IopScfKeys.hh>
#include <molsturm/RestrictedClosedIntegralOperator.hh>
#include <molsturm/ScfDebugWrapper.hh>
#include <molsturm/scf_guess.hh>
#include <molsturm/version.hh>

DefException1(
      ExcTooSmallBasis, size_t,
      << "A basis of size " << arg1
      << " is too small to incorporate this many electrons. Choose a larger basis.");

namespace hf {
using namespace molsturm;
using namespace linalgwrap;
using namespace krims;

/** Print the resulting eigenstates */
template <typename State>
void print_res(const State& res) {
  io::OstreamState oldstate(std::cout);

  const std::string indt = "    ";
  std::cout << "Obtained Hartree-Fock orbitals: " << std::endl;
  std::cout << indt << "a b" << std::endl;

  size_t n_alpha = res.problem_matrix().n_alpha();
  size_t n_beta = res.problem_matrix().n_beta();
  const auto& orben = res.eigensolution().evalues();
  for (const auto& val : orben) {
    std::cout << indt;
    if (n_alpha > 0) {
      --n_alpha;
      std::cout << "* ";
    } else {
      std::cout << "  ";
    }
    if (n_beta > 0) {
      --n_beta;
      std::cout << "* ";
    } else {
      std::cout << "  ";
    }
    std::cout << std::right << std::setw(17) << std::setprecision(12) << val << std::endl;
  }
}

/** Run an (atomic) SCF based on the integral data Sturmian14.
 *
 * \param Z       Number of nuclei
 * \param k_exp   Sturmian exponent
 * \param n_alpha Number of alpha electrons
 * \param n_beta  Number of beta electrons
 */
void run_rhf_sturmian(args_type args, bool debug = false) {
  //
  // Types and settings
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;

  // The lookup class type to get the actual integrals
  typedef gint::IntegralLookup<gint::COMPLEX_ATOMIC> int_lookup_type;

  // The type of the integral terms:
  typedef typename int_lookup_type::integral_type integral_type;

  //
  // Integral terms
  //
  // Generate integral lookup object
  krims::ParameterMap intparams{
        {"basis_type", args.basis_type},
        {"k_exponent", args.k_exp},
        {"Z_charge", args.Z_charge},
        {"n_max", static_cast<int>(args.n_max)},
        {"l_max", static_cast<int>(args.n_max) - 1},
        {"m_max", static_cast<int>(args.n_max) - 1},
  };
  int_lookup_type integrals{intparams};

  // Get the integral as actual objects.
  integral_type S_bb = integrals("overlap");
  integral_type T_bb = integrals("kinetic");
  integral_type V0_bb = integrals("nuclear_attraction");
  integral_type J_bb = integrals("coulomb");
  integral_type K_bb = integrals("exchange");

  std::cout << "Basis size:  " << S_bb.n_rows() << std::endl << std::endl;

  // Combine 1e terms:
  std::vector<integral_type> terms_1e{std::move(T_bb), std::move(V0_bb)};

  //
  // Problem setup
  //
  assert_throw(args.n_eigenpairs >= std::max(args.n_alpha, args.n_beta),
               krims::ExcTooLarge<size_t>(std::max(args.n_alpha, args.n_beta),
                                          args.n_eigenpairs));
  assert_throw(2 * std::max(args.n_alpha, args.n_beta) + 2 <= S_bb.n_rows(),
               ExcTooSmallBasis(S_bb.n_rows()));

  // The term container for the fock operator matrix
  IntegralTermContainer<stored_matrix_type> integral_container(
        std::move(terms_1e), std::move(J_bb), std::move(K_bb));

  RestrictedClosedIntegralOperator<stored_matrix_type> fock_bb(integral_container,
                                                               args.n_alpha, args.n_beta);

  // Update with a guess solution
  auto guess_solution =
        scf_guess(fock_bb, S_bb, {{ScfGuessKeys::method, /*"hcore"*/ "loewdin"}});
  fock_bb.update(guess_solution.evectors_ptr);

  krims::ParameterMap params{
        // error
        {IopScfKeys::max_error_norm, args.error},
        {IopScfKeys::max_1e_energy_change, args.error * 100.},
        {IopScfKeys::max_tot_energy_change, args.error / 4.},
        //
        {IopScfKeys::max_iter, 25ul},
        {IopScfKeys::n_eigenpairs, args.n_eigenpairs},
        {IopScfKeys::verbosity, ScfMsgType::FinalSummary | ScfMsgType::IterationProcess},
        {IopScfKeys::n_prev_steps, args.diis_size},
  };

  if (debug) {
    std::ofstream mathematicafile("/tmp/debug_molsturm_rhf_sturmian.m");
    auto debugout =
          linalgwrap::io::make_formatted_stream_writer<linalgwrap::io::Mathematica,
                                                       scalar_type>(mathematicafile,
                                                                    1e-12);
    debugout.write("guess", guess_solution.evectors());
    debugout.write("sbb", S_bb);
    typedef gscf::PulayDiisScfState<decltype(fock_bb), decltype(S_bb)> inner_state_type;
    typedef molsturm::detail::IopScfStateWrapper<inner_state_type> state_type;
    molsturm::detail::IopScfWrapper<gscf::PulayDiisScf<state_type>> solver{params};
    ScfDebugWrapper<decltype(solver)> solwrap(debugout, solver);
    auto res = solwrap.solve(fock_bb, S_bb);
    print_res(res);
  } else {
    auto res = run_scf(fock_bb, S_bb, params, guess_solution);
    print_res(res);
  }
}

int main(int argc, char** argv) {
  std::cout << "molsturm version:   " << molsturm::version::version_string() << std::endl
            << "gscf version:       " << gscf::version::version_string() << std::endl
            << "gint version:       " << gint::version::version_string() << std::endl
            << "linalgwrap version: " << linalgwrap::version::version_string()
            << std::endl;

  args_type args;
  if (!parse_args(argc, argv, args)) return 1;
  std::cout << "The following configuration was read:" << std::endl << args << std::endl;

  run_rhf_sturmian(args, /*debug = */ false);
  return 0;
}

}  // namespace hf

int main(int argc, char** argv) { return hf::main(argc, argv); }
