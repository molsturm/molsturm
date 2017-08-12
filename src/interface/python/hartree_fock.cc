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

#include "hartree_fock.hh"
#include "ExcInvalidParameters.hh"
#include "config.hh"
#include "export_hf_results.hh"
#include <gint/IntegralLookup.hh>
#include <gint/IntegralType.hh>
#include <molsturm/FockOperator.hh>
#include <molsturm/IopScf.hh>
#include <molsturm/MolecularSystem.hh>
#include <molsturm/OverlapMatrix.hh>
#include <molsturm/scf_guess.hh>

namespace molsturm {
namespace iface {

template <RestrictionType restrict>
HfResults hartree_fock_inner(const HfParameters& params) {
  // Lookup integral terms
  gint::IntegralLookup<matrix_type> integrals(params.integral_params);
  auto Sa_bb = integrals.lookup_integral(gint::IntegralTypeKeys::overlap);

  // Checks about basis size:
  const size_t n_bas        = Sa_bb.n_rows();
  const size_t max_elec     = std::max(params.n_alpha, params.n_beta);
  const size_t n_eigenpairs = params.scf_params.at<size_t>("n_eigenpairs");
  assert_throw(
        max_elec < n_bas,
        ExcInvalidParameters("A basis with n_bas == " + std::to_string(n_bas) +
                             " basis functions is too small for treating a system with " +
                             std::to_string(params.n_alpha) + " alpha and " +
                             std::to_string(params.n_beta) + " beta electrons."));
  assert_throw(
        n_eigenpairs >= 2 * max_elec,
        ExcInvalidParameters(
              "Cannot treat a system with " + std::to_string(params.n_alpha) +
              " alpha and " + std::to_string(params.n_beta) +
              " beta electrons with computing only " + std::to_string(n_eigenpairs) +
              " eigenpairs in the SCF eigensolver. You need to request at least " +
              std::to_string(2 * max_elec) + " eigenpairs."));

  // Run solver
  OverlapMatrix<matrix_type, restrict> S_bb(Sa_bb);
  MolecularSystem system(params.structure, {{params.n_alpha, params.n_beta}});
  FockOperator<matrix_type, restrict> fock_bb(integrals, system);
  auto guess = scf_guess(system, fock_bb, S_bb, params.guess_params);
  fock_bb.update(guess.evectors_ptr);

  auto result = run_scf(fock_bb, S_bb, guess, params.scf_params);
  return export_hf_results(result, integrals.eri_tensor(), params);
}

HfResults restricted_hartree_fock(const HfParameters& p) {
  return hartree_fock_inner<RestrictionType::RestrictedClosed>(p);
}

HfResults unrestricted_hartree_fock(const HfParameters& p) {
  return hartree_fock_inner<RestrictionType::Unrestricted>(p);
}

}  // namespace iface
}  // namespace molsturm
