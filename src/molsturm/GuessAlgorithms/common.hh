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
#include "GuessAlgorithmsKeysBase.hh"
#include <gint/IntegralType.hh>
#include <gscf/FocklikeMatrix_i.hh>
#include <lazyten/eigensystem.hh>
#include <lazyten/rescue.hh>
#include <molsturm/MolecularSystem.hh>

namespace molsturm {

/** Thrown if obtaining the scf guess failed (mostly because the eigensolver
 * had trouble finding the required eigenvalues)
 */
DefException1(ExcObtainingScfGuessFailed, std::string,
              << "Obtaining the guess for the SCF failed. Reason: " << arg1);

/** Thrown if solver parameters are wrong or inconsistent */
DefException1(ExcInvalidScfGuessParametersEncountered, std::string,
              << "The Scf guess could not be obtained, because the set of "
                 "parameters passed is not valid. Details: "
              << arg1);

/** Take separate guesses for the alpha and the beta block
 * and combine them to a single solution by padding with
 * zeros.
 *
 * This will make the guess suitable for unrestricted calculations
 */
template <typename Solution>
Solution make_block_diagonal(const Solution& solution_alpha,
                             const Solution& solution_beta) {
  typedef typename Solution::evector_type evector_type;

  Solution solution;
  auto& evecs_alpha   = solution_alpha.evectors();
  auto& evecs_beta    = solution_beta.evectors();
  const size_t n_orbs = evecs_alpha.n_vectors();
  const size_t n_bas  = evecs_alpha.n_elem();
  assert_size(n_orbs, evecs_beta.n_vectors());
  assert_size(n_bas, evecs_beta.n_elem());

  lazyten::MultiVector<evector_type> guess(2 * n_bas, 2 * n_orbs);
  for (size_t f = 0; f < n_orbs; ++f) {
    std::copy(evecs_alpha[f].begin(), evecs_alpha[f].end(), guess[f].begin());
    std::copy(evecs_beta[f].begin(), evecs_beta[f].end(),
              guess[f + n_orbs].begin() + n_bas);
  }
  solution.evectors() = guess;

  auto& evals_alpha = solution_alpha.evalues();
  auto& evals_beta  = solution_beta.evalues();
  solution.evalues().resize(2 * n_orbs);
  auto it = std::copy(evals_alpha.begin(), evals_alpha.end(), solution.evalues().begin());
  it      = std::copy(evals_beta.begin(), evals_beta.end(), it);
  assert_internal(it == solution.evalues().end());

  return solution;
}

/** Replicate a guess solution for a restricted closed operator
 * such that it exists twice on the diagonal.
 * with zeros padded before and after in the coefficients.
 *
 * This will make the guess suitable for unrestricted calculations
 */
template <typename Solution>
Solution replicate_block(const Solution& solution) {
  return make_block_diagonal(solution, solution);
}

}  // namespace molsturm
