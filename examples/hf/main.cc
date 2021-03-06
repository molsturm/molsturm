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

#define KRIMS_INIT_EXCEPTION_SYSTEM
#include <krims/ExceptionSystem.hh>

#include "parse_args.hh"

#include <gint/IntegralLookup.hh>
#include <gint/IntegralLookupKeys.hh>
#include <gint/OrbitalType.hh>
#include <gint/version.hh>
#include <gscf/version.hh>
#include <iostream>
#include <lazyten/SmallMatrix.hh>
#include <lazyten/io.hh>
#include <lazyten/version.hh>
#include <molsturm/FockOperator.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/IopScfKeys.hh>
#include <molsturm/OverlapMatrix.hh>
#include <molsturm/ScfDebugWrapper.hh>
#include <molsturm/Version.hh>
#include <molsturm/scf_guess.hh>

DefException1(
      ExcTooSmallBasis, size_t,
      << "A basis of size " << arg1
      << " is too small to incorporate this many electrons. Choose a larger basis.");

namespace hf {
using namespace molsturm;
using namespace lazyten;
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
template <RestrictionType restrict>
void run_hf(args_type args, bool debug = false) {
  using gint::IntegralLookupKeys;
  using gint::IntegralTypeKeys;

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
          {"basis_set_name", args.basis_set},
          //
          {IntegralLookupKeys::structure, args.system.structure},
    });
  } else {
    assert_implemented(false);
  }

  integral_lookup_type integrals(std::move(intparams));
  auto Sa_bb = integrals.lookup_integral(IntegralTypeKeys::overlap);

  //
  // Checks about basis size:
  //
  std::cout << "Basis size:  " << Sa_bb.n_rows() << std::endl << std::endl;

  const size_t max_elec = std::max(args.system.n_alpha, args.system.n_beta);
  assert_throw(max_elec < Sa_bb.n_rows(), ExcTooSmallBasis(Sa_bb.n_rows()));

  assert_throw(2 * args.n_eigenpairs >= std::max(args.system.n_alpha, args.system.n_beta),
               krims::ExcTooLarge<size_t>(max_elec, 2 * args.n_eigenpairs));

  //
  // Problem setup
  //
  OverlapMatrix<stored_matrix_type, restrict> S_bb(Sa_bb);
  FockOperator<stored_matrix_type, restrict> fock_bb(integrals, args.system);

  // Obtain an SCF guess
  krims::GenMap guess_params{{ScfGuessKeys::method, args.guess_method}};
  guess_params.update(ScfGuessKeys::eigensolver_params,
                      {{EigensystemSolverKeys::method, args.guess_esolver}});
  auto guess = scf_guess(args.system, fock_bb, S_bb, guess_params);

  const size_t n_eigenpairs = restrict == RestrictionType::RestrictedClosed
                                    ? args.n_eigenpairs / 2
                                    : args.n_eigenpairs;
  krims::GenMap params{
        // error
        {IopScfKeys::max_error_norm, args.error},
        {IopScfKeys::max_1e_energy_change, args.error * 100.},
        {IopScfKeys::max_tot_energy_change, args.error / 4.},
        //
        {IopScfKeys::max_iter, args.max_iter},
        {IopScfKeys::n_eigenpairs, n_eigenpairs},
        {IopScfKeys::verbosity, ScfMsgType::FinalSummary | ScfMsgType::IterationProcess},
        {gscf::PulayDiisScfKeys::n_prev_steps, args.diis_size},
  };
  params.update(IopScfKeys::eigensolver_params,
                {{EigensystemSolverKeys::method, args.eigensolver}});

  if (debug) {
    std::ofstream mathematicafile("/tmp/debug_molsturm_rhf_sturmian.m");
    auto debugout = lazyten::io::make_writer<lazyten::io::Mathematica, scalar_type>(
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
  std::cout << "molsturm version: " << molsturm::Version::as_string() << std::endl
            << "gscf version:     " << gscf::version::version_string() << std::endl
            << "gint version:     " << gint::version::version_string() << std::endl
            << "lazyten version:  " << lazyten::version::version_string() << std::endl;

  args_type args;
  if (!parse_args(argc, argv, args)) return 1;

  const int prec = static_cast<int>(std::cout.precision());
  std::cout << '\n'
            << "The following configuration was read:\n"
            << std::setprecision(15) << args << std::setprecision(prec) << std::endl;

  const bool debug = false;
  if (args.restricted) {
    run_hf<RestrictionType::RestrictedClosed>(args, debug);
  } else {
    run_hf<RestrictionType::Unrestricted>(args, debug);
  }
  return 0;
}

}  // namespace hf

int main(int argc, char** argv) { return hf::main(argc, argv); }
