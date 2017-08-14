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

  // Build the molecular system
  const krims::GenMap system_params = params.submap("system");
  const MolecularSystem system{system_params.at<gint::Structure>(SystemKeys::structure),
                               {{system_params.at<size_t>(SystemKeys::n_alpha),
                                 system_params.at<size_t>(SystemKeys::n_beta)}}};

  // Overlap matrix and  Fock or Kohn-Sham matrix setup
  auto Sa_bb = integrals.lookup_integral(gint::IntegralTypeKeys::overlap);
  OverlapMatrix<matrix_type, restrict> S_bb(Sa_bb);
  ScfOperator scfop_bb(integrals, system);

  // Obtaining and setting the guess.
  const krims::GenMap guess_params = params.submap("guess");
  const std::string guess_method   = guess_params.at<std::string>(ScfGuessKeys::method);

  auto guess = (guess_method == std::string("external"))
                     ? solution_view.as_eigensolution()
                     : scf_guess(system, scfop_bb, S_bb, guess_params);
  scfop_bb.update(guess.evectors_ptr);

  // Run SCF and export results back to caller.
  auto result = run_scf(scfop_bb, S_bb, guess, params.submap("scf"));
  solution_view.set_from_eigensolution(result.eigensolution());
  return export_hf_results(result, integrals.eri_tensor());
}

ScfResults self_consistent_field(const ScfKind type, const ScfParameters& params,
                                 ScfSolutionView& solution_view) {
  const bool restricted = type == ScfKind::RHF;

  // Add the structure for the integral library if not done yet:
  params.insert_default("integrals/structure", params.at_raw_value("system/structure"));

  // TODO Potentially setting up this thing is expensive.
  //      Ideally we want to keep this around until it is no longer needed
  //      in other words, we probably want to pass it back to the caller
  //      inside the ScfResults or so.
  gint::IntegralLookup<matrix_type> integrals(params.submap("integrals"));
  const size_t n_bas =
        integrals.lookup_integral(gint::IntegralTypeKeys::overlap).n_rows();
  const size_t n_alpha  = params.at<size_t>("system/" + SystemKeys::n_alpha);
  const size_t n_beta   = params.at<size_t>("system/" + SystemKeys::n_beta);
  const size_t max_elec = std::max(n_alpha, n_beta);
  const size_t n_spin   = solution_view.n_spin;

  //
  // Size checks
  //
  assert_throw((restricted && n_spin == 1) || (!restricted && n_spin == 2),
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
                             std::to_string(n_alpha) + " alpha and " +
                             std::to_string(n_beta) + " beta electrons."));

  assert_throw(max_elec < solution_view.n_fock,
               ExcInvalidParameters("Requested too few eigeinpairs. Cannot incorporate " +
                                    std::to_string(n_alpha) + " alpha and " +
                                    std::to_string(n_beta) + " beta electrons if only " +
                                    std::to_string(solution_view.n_fock) +
                                    " orbitals are computed for each spin."));

  //
  // Make sure the n_eigenpairs setting exists and the number of eigenpairs we
  // can store in the arrays agrees with it.
  //

  // Note: The insert_default only does something if the key does *not* exist inside the
  //       params yet.
  const std::string eigenpairs_key = "scf/" + gscf::ScfBaseKeys::n_eigenpairs;
  params.insert_default(eigenpairs_key, solution_view.n_fock * solution_view.n_spin);
  const size_t n_eigenpairs = params.at<size_t>(eigenpairs_key);

  assert_throw(n_eigenpairs == solution_view.n_fock * n_spin,
               ExcInvalidParameters(
                     "Number of eigenpairs requested via the ScfParameters (== " +
                     std::to_string(n_eigenpairs) +
                     ") and the size of the arrays "
                     "which is reserved for the "
                     "energies/coeffients in the ScfSolutionView (== " +
                     std::to_string(solution_view.n_fock) + " orbitals * " +
                     std::to_string(n_spin) + " spin components) does not agree."));

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
