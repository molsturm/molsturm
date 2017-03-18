#define KRIMS_INIT_EXCEPTION_SYSTEM
#include <krims/ExceptionSystem.hh>

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

  const std::string ind = "    ";
  std::cout << "Obtained Hartree-Fock orbitals: " << std::endl;
  std::cout << ind << "a b" << std::endl;

  auto occa = res.problem_matrix().indices_subspace(gscf::OrbitalSpace::OCC_ALPHA);
  auto occb = res.problem_matrix().indices_subspace(gscf::OrbitalSpace::OCC_BETA);
  assert_dbg(occa == occb, krims::ExcNotImplemented());

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

  //
  // Types and settings
  //
  // Types of scalar and matrix
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef gint::Integral<stored_matrix_type> integral_type;

  // The lookup class type to get the actual integrals
  // TODO This whole handling would be much easier if the orbital type was not in there!
  typedef gint::IntegralLookup<gint::OrbitalType::COMPLEX_ATOMIC> int_lookup_ca_type;
  typedef gint::IntegralLookup<gint::OrbitalType::REAL_MOLECULAR> int_lookup_rm_type;
  std::unique_ptr<int_lookup_ca_type> integrals_ca_ptr;
  std::unique_ptr<int_lookup_rm_type> integrals_rm_ptr;

  //
  // Lookup integral terms
  //
  // The integral objects to be filled:
  integral_type S_bb, T_bb, V0_bb, J_bb, K_bb;

  if (args.sturmian) {
    assert_throw(args.system.structure.n_atoms() == 1, krims::ExcNotImplemented());

    krims::GenMap intparams{
          {"basis_type", args.basis_type},
          {"Z_charge", args.system.structure.total_charge()},
          {"structure", args.system.structure},
          {"k_exponent", args.k_exp},
          {"n_max", static_cast<int>(args.n_max)},
          {"l_max", static_cast<int>(args.l_max)},
          {"m_max", static_cast<int>(args.m_max)},
    };
    integrals_ca_ptr.reset(new int_lookup_ca_type(std::move(intparams)));

    S_bb = integrals_ca_ptr->lookup_integral(IntegralTypeKeys::overlap);
    T_bb = integrals_ca_ptr->lookup_integral(IntegralTypeKeys::kinetic);
    V0_bb = integrals_ca_ptr->lookup_integral(IntegralTypeKeys::nuclear_attraction);
    J_bb = integrals_ca_ptr->lookup_integral(IntegralTypeKeys::coulomb);
    K_bb = integrals_ca_ptr->lookup_integral(IntegralTypeKeys::exchange);

    const size_t max_elec = std::max(args.system.n_alpha, args.system.n_beta);
    assert_throw(2 * max_elec + 2 <= S_bb.n_rows(), ExcTooSmallBasis(S_bb.n_rows()));
  } else if (args.gaussian) {
    krims::GenMap intparams{
          {"basis_type", args.basis_type},
          {"structure", args.system.structure},
          {"basis_set", args.basis_set},
    };
    integrals_rm_ptr.reset(new int_lookup_rm_type(std::move(intparams)));

    S_bb = integrals_rm_ptr->lookup_integral(IntegralTypeKeys::overlap);
    T_bb = integrals_rm_ptr->lookup_integral(IntegralTypeKeys::kinetic);
    V0_bb = integrals_rm_ptr->lookup_integral(IntegralTypeKeys::nuclear_attraction);
    J_bb = integrals_rm_ptr->lookup_integral(IntegralTypeKeys::coulomb);
    K_bb = integrals_rm_ptr->lookup_integral(IntegralTypeKeys::exchange);
  }

  //
  // Print basis info:
  //
  std::cout << "Basis size:  " << S_bb.n_rows() << std::endl << std::endl;
  assert_throw(
        args.n_eigenpairs >= std::max(args.system.n_alpha, args.system.n_beta),
        krims::ExcTooLarge<size_t>(std::max(args.system.n_alpha, args.system.n_beta),
                                   args.n_eigenpairs));

  //
  // Problem setup
  //
  // The term container for the fock operator matrix
  IntegralTermContainer<stored_matrix_type> integral_container(
        {{std::move(T_bb), std::move(V0_bb)}}, std::move(J_bb), std::move(K_bb));

  RestrictedClosedIntegralOperator<stored_matrix_type> fock_bb(integral_container,
                                                               args.system);

  // Obtain an SCF guess
  auto guess = scf_guess(args.system, fock_bb, S_bb,
                         {{ScfGuessKeys::method, args.guess_method}});

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
  std::cout << "molsturm version:   " << molsturm::version::version_string() << std::endl
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
