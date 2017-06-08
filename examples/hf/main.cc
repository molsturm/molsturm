#define KRIMS_INIT_EXCEPTION_SYSTEM
#include <krims/ExceptionSystem.hh>

#include "parse_args.hh"

#include <gint/IntegralLookup.hh>
#include <gint/IntegralLookupKeys.hh>
#include <gint/OrbitalType.hh>
#include <gint/version.hh>
#include <gscf/version.hh>
#include <iostream>
#include <linalgwrap/io.hh>
#include <linalgwrap/version.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/IopScfKeys.hh>
#include <molsturm/RestrictedClosedIntegralOperator.hh>
#include <molsturm/ScfDebugWrapper.hh>
#include <molsturm/Version.hh>
#include <molsturm/scf_guess.hh>

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

  const std::string ind = "    ";
  std::cout << "Obtained Hartree-Fock orbitals: " << std::endl;
  std::cout << ind << "a b" << std::endl;

  auto occa = res.problem_matrix().indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
  auto occb = res.problem_matrix().indices_orbspace(gscf::OrbitalSpace::OCC_BETA);
  assert_implemented(occa == occb);

  const auto& orben = res.eigensolution().evalues();
  for (size_t i = 0; i < orben.size(); ++i) {
    std::cout << ind;
    std::cout << (occa.contains(i) ? '*' : ' ');
    std::cout << ' ';
    std::cout << (occb.contains(i) ? '*' : ' ');
    std::cout << ' ';
    std::cout << std::right << std::setw(17) << std::setprecision(12) << orben[i]
              << std::endl;
  }
}

/** Run an SCF */
void run_rhf(args_type args, bool debug = false) {
  using gint::IntegralTypeKeys;
  using gint::IntegralLookupKeys;

  //
  // Types and settings
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef gint::IntegralLookup<stored_matrix_type> integral_lookup_type;

  //
  // Lookup integral terms
  //
  krims::GenMap intparams;

  if (args.sturmian) {
    assert_implemented(args.system.structure.n_atoms() == 1);
    intparams.update({
          {IntegralLookupKeys::basis_type, args.basis_type},
          {IntegralLookupKeys::orbital_type, gint::OrbitalType::COMPLEX_ATOMIC},
          //
          {"k_exponent", args.k_exp},
          //
          {IntegralLookupKeys::structure, args.system.structure},
    });

    if (args.nlm_basis.size() > 0) {
      intparams.update("nlm_basis", args.nlm_basis);
    } else {
      intparams.update({
            {"n_max", static_cast<int>(args.n_max)},
            {"l_max", static_cast<int>(args.l_max)},
            {"m_max", static_cast<int>(args.m_max)},
      });
    }
  } else if (args.gaussian) {
    intparams.update({
          {IntegralLookupKeys::basis_type, args.basis_type},
          {IntegralLookupKeys::orbital_type, gint::OrbitalType::REAL_MOLECULAR},
          {"basis_set", args.basis_set},
          //
          {IntegralLookupKeys::structure, args.system.structure},
    });
  } else {
    assert_implemented(false);
  }

  integral_lookup_type integrals(std::move(intparams));
  auto S_bb = integrals.lookup_integral(IntegralTypeKeys::overlap);
  auto T_bb = integrals.lookup_integral(IntegralTypeKeys::kinetic);
  auto V0_bb = integrals.lookup_integral(IntegralTypeKeys::nuclear_attraction);
  auto J_bb = integrals.lookup_integral(IntegralTypeKeys::coulomb);
  auto K_bb = integrals.lookup_integral(IntegralTypeKeys::exchange);

  //
  // Checks about basis size:
  //
  std::cout << "Basis size:  " << S_bb.n_rows() << std::endl << std::endl;

  const size_t max_elec = std::max(args.system.n_alpha, args.system.n_beta);
  assert_throw(max_elec < S_bb.n_rows(), ExcTooSmallBasis(S_bb.n_rows()));

  assert_throw(args.n_eigenpairs >= std::max(args.system.n_alpha, args.system.n_beta),
               krims::ExcTooLarge<size_t>(max_elec, args.n_eigenpairs));

  //
  // Problem setup
  //
  // The term container for the fock operator matrix
  IntegralTermContainer<stored_matrix_type> integral_container(
        {{std::move(T_bb), std::move(V0_bb)}}, std::move(J_bb), std::move(K_bb));

  RestrictedClosedIntegralOperator<stored_matrix_type> fock_bb(integral_container,
                                                               args.system);

  // Obtain an SCF guess
  krims::GenMap guess_params{{ScfGuessKeys::method, args.guess_method}};
  guess_params.update(ScfGuessKeys::eigensolver_params,
                      {{EigensystemSolverKeys::method, args.guess_esolver}});
  auto guess = scf_guess(args.system, fock_bb, S_bb, guess_params);

  krims::GenMap params{
        // error
        {IopScfKeys::max_error_norm, args.error},
        {IopScfKeys::max_1e_energy_change, args.error * 100.},
        {IopScfKeys::max_tot_energy_change, args.error / 4.},
        //
        {IopScfKeys::max_iter, args.max_iter},
        {IopScfKeys::n_eigenpairs, args.n_eigenpairs},
        {IopScfKeys::verbosity, ScfMsgType::FinalSummary | ScfMsgType::IterationProcess},
        {gscf::PulayDiisScfKeys::n_prev_steps, args.diis_size},
  };
  params.update(IopScfKeys::eigensolver_params,
                {{EigensystemSolverKeys::method, args.eigensolver}});

  if (debug) {
    std::ofstream mathematicafile("/tmp/debug_molsturm_rhf_sturmian.m");
    auto debugout = linalgwrap::io::make_writer<linalgwrap::io::Mathematica, scalar_type>(
          mathematicafile, 1e-20);
    debugout.write("guess", guess.evectors());
    debugout.write("sbb", S_bb);
    typedef gscf::PulayDiisScfState<decltype(fock_bb), decltype(S_bb)> inner_state_type;
    typedef molsturm::detail::IopScfStateWrapper<inner_state_type> state_type;
    molsturm::detail::IopScfWrapper<gscf::PulayDiisScf<state_type>> solver{params};
    ScfDebugWrapper<decltype(solver)> solwrap(debugout, solver);

    fock_bb.update(guess.evectors_ptr);
    auto res = solwrap.solve(fock_bb, S_bb);
    print_res(res);
  } else {
    fock_bb.update(guess.evectors_ptr);
    auto res = run_scf(fock_bb, S_bb, guess, params);
    print_res(res);
  }
}

int main(int argc, char** argv) {
  std::cout << "molsturm version:   " << molsturm::Version::as_string() << std::endl
            << "gscf version:       " << gscf::version::version_string() << std::endl
            << "gint version:       " << gint::version::version_string() << std::endl
            << "linalgwrap version: " << linalgwrap::version::version_string()
            << std::endl;

  args_type args;
  if (!parse_args(argc, argv, args)) return 1;

  const int prec = static_cast<int>(std::cout.precision());
  std::cout << '\n'
            << "The following configuration was read:\n"
            << std::setprecision(15) << args << std::setprecision(prec) << std::endl;

  run_rhf(args, /*debug = */ false);
  return 0;
}

}  // namespace hf

int main(int argc, char** argv) { return hf::main(argc, argv); }
