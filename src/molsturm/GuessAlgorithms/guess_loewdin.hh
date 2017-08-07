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

#pragma once
#include "common.hh"

namespace molsturm {

typedef GuessAlgorithmsKeysBase GuessLoewdinKeys;

/** Obtain a Löwdin guess: Diagonalise the overlap
 *
 * If obtaining the guess failed, the function throws
 * ExcObtainingScfGuessFailed.
 *
 * ## Parameters
 * The algorithm takes some parameters, which change the internal behaviour.
 * For the full list see the class GuessLoewdinKeys
 *
 * - "eigensolver":   Parameter submap which controls the behaviour of
 *                    the eigensolver. By default the solver will try to
 *                    find the smallest real eigenvalues
 *
 * \throws ExcObtainingScfGuessFailed if obtaining the guess failed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
lazyten::EigensolutionTypeFor<true, IntegralOperator> guess_loewdin(
      const MolecularSystem& /*system*/, const IntegralOperator& fock_bb,
      const OverlapMatrix& S_bb, const krims::GenMap& params) {
  using namespace lazyten;

  // apply Löwdin normalisation to the basis functions
  //   - Diagonalise the overlap
  //   - Take 1/\sqrt{evals} at the diagonal
  //   - results in orthonormalised basis functions

  krims::GenMap eigensolver_params = params.submap(GuessLoewdinKeys::eigensolver_params);
  eigensolver_params.insert_default(EigensystemSolverKeys::which, "LR");

  const auto occa        = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
  const auto occb        = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA);
  const size_t n_vectors = std::max(occa.size(), occb.size());

  // Get alpha-alpha block of the overlap matrix.
  // Note by construction this block is identical to the beta-beta block
  const auto& Sa_bb = S_bb.block_alpha();

  // Restricted open-shell is not yet implemented
  assert_internal(occa == occb || !fock_bb.restricted());
  try {
    auto sol = eigensystem_hermitian(Sa_bb, n_vectors, eigensolver_params);

    // Eigenvectors and eigenvalues.
    auto& evectors = sol.evectors();
    auto& evalues  = sol.evalues();

    assert_internal(evectors.n_vectors() == n_vectors);
    assert_internal(evectors.n_elem() == Sa_bb.n_cols());

    for (size_t i = 0; i < evectors.n_vectors(); ++i) {
      evectors[i] *= 1. / sqrt(evalues[i]);
      evalues[i] = 1.;
    }
    return fock_bb.restricted() ? sol : replicate_block(sol);
  } catch (const SolverException& e) {
    rescue::failed_eigenproblem(
          Eigenproblem<true, typename OverlapMatrix::int_term_type>(Sa_bb), params);
    assert_throw(false,
                 ExcObtainingScfGuessFailed(
                       "Eigensolver for overlap failed with message " + e.extra()));
    return EigensolutionTypeFor<true, IntegralOperator>{};
  }
}

}  // namespace molsturm
