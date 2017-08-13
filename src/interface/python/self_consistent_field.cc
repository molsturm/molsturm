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

#include "self_consistent_field.hh"
#include "ExcInvalidParameters.hh"
#include "ScfSolutionView.hh"
#include "config.hh"
#include <gint/IntegralLookup.hh>
#include <krims/ExceptionSystem.hh>
#include <molsturm/FockOperator.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/MolecularSystem.hh>
#include <molsturm/OverlapMatrix.hh>
#include <molsturm/scf_guess.hh>

#include "export_hf_results.hh"

namespace molsturm {
namespace iface {

template <typename ScfOperator>
ScfResults scf_for_operator(const ScfParameters& params,
                            const gint::IntegralLookup<matrix_type>& integrals,
                            ScfSolutionView& solution_view) {
  constexpr RestrictionType restrict = ScfOperator::restriction_type;
  MolecularSystem system(params.structure, {{params.n_alpha, params.n_beta}});

  // Overlap matrix and  Fock or Kohn-Sham matrix setup
  auto Sa_bb = integrals.lookup_integral(gint::IntegralTypeKeys::overlap);
  OverlapMatrix<matrix_type, restrict> S_bb(Sa_bb);
  ScfOperator scfop_bb(integrals, system);

  // Obtaining and setting the guess.
  const std::string guess_method =
        params.scf_params.at<std::string>(ScfGuessKeys::method);
  auto guess = (guess_method == std::string("external"))
                     ? solution_view.as_eigensolution()
                     : scf_guess(system, scfop_bb, S_bb, params.guess_params);
  scfop_bb.update(guess.evectors_ptr);

  // Run SCF and export results back to caller.
  auto result = run_scf(scfop_bb, S_bb, guess, params.scf_params);
  solution_view.set_from_eigensolution(result.eigensolution());
  return export_hf_results(result, integrals.eri_tensor(), params);
}

ScfResults self_consistent_field(const ScfKind type, const ScfParameters& params,
                                 ScfSolutionView& solution_view) {
  const bool restricted = type == ScfKind::RHF;

  // TODO Potentially setting up this thing is expensive.
  //      Ideally we want to keep this around until it is no longer needed
  //      in other words, we probably want to pass it back to the caller
  //      inside the ScfResults or so.
  gint::IntegralLookup<matrix_type> integrals(params.integral_params);
  const size_t n_bas =
        integrals.lookup_integral(gint::IntegralTypeKeys::overlap).n_rows();
  const size_t max_elec = std::max(params.n_alpha, params.n_beta);

  //
  // Size checks
  //
  assert_throw((restricted && solution_view.n_spin == 1) ||
                     (!restricted && solution_view.n_spin == 2),
               ExcInvalidParameters("For restricted SCF we need exactly 1 spin component "
                                    "and for unrestricted SCF exactly 2."));

  assert_throw(n_bas == solution_view.n_bas,
               ExcInvalidParameters("The number of basis functions inferred via the "
                                    "orbcoeff_bf array of the ScfSolutionView (== " +
                                    std::to_string(solution_view.n_bas) +
                                    ") and the number of basis functions which are "
                                    "supplied via the integral backend (== " +
                                    std::to_string(n_bas) + ") differs."));

  assert_throw(
        max_elec < n_bas,
        ExcInvalidParameters("A basis with n_bas == " + std::to_string(n_bas) +
                             " basis functions is too small for treating a system with " +
                             std::to_string(params.n_alpha) + " alpha and " +
                             std::to_string(params.n_beta) + " beta electrons."));

  // Make sure the n_eigenpairs setting and the number of eigenpairs we can store in the
  // arrays agrees.
  if (params.scf_params.exists(gscf::ScfBaseKeys::n_eigenpairs)) {
    const size_t n_eigenpairs =
          params.scf_params.at<size_t>(gscf::ScfBaseKeys::n_eigenpairs);
    assert_throw(n_eigenpairs == solution_view.n_fock * solution_view.n_spin,
                 ExcInvalidParameters(
                       "Number of eigenpairs requested via the ScfParameters (== " +
                       std::to_string(n_eigenpairs) +
                       ") and the size of the arrays "
                       "which is reserved for the "
                       "energies/coeffients in the ScfSolutionView (== " +
                       std::to_string(solution_view.n_fock) + " orbitals * " +
                       std::to_string(solution_view.n_spin) +
                       " spin components) does not agree."));
  } else {
    params.scf_params.insert_default(gscf::ScfBaseKeys::n_eigenpairs,
                                     solution_view.n_fock * solution_view.n_spin);
  }

  const size_t n_eigenpairs =
        params.scf_params.at<size_t>(gscf::ScfBaseKeys::n_eigenpairs);
  assert_throw(
        n_eigenpairs >= 2 * max_elec,
        ExcInvalidParameters(
              "Cannot treat a system with " + std::to_string(params.n_alpha) +
              " alpha and " + std::to_string(params.n_beta) +
              " beta electrons with computing only " + std::to_string(n_eigenpairs) +
              " eigenpairs in the SCF eigensolver. You need to request at least " +
              std::to_string(2 * max_elec) + " eigenpairs."));

  //
  // Run the SCF
  //
  typedef FockOperator<matrix_type, RestrictionType::RestrictedClosed> fock_rhf_t;
  typedef FockOperator<matrix_type, RestrictionType::Unrestricted> fock_uhf_t;

  switch (type) {
    case ScfKind::RHF:
      return scf_for_operator<fock_rhf_t>(params, integrals, solution_view);
    case ScfKind::UHF:
      return scf_for_operator<fock_uhf_t>(params, integrals, solution_view);
    default:
      assert_throw(false,
                   ExcInvalidParameters(
                         "The ScfKind / SCF type parameter supplied is not understood."));
  }
}

}  // namespace iface
}  // namespace molsturm
