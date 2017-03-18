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

typedef GuessAlgorithmsKeysBase GuessFromArchivedKeys;

/** Obtain a guess by extrapolating based on an archived solution
 *
 * ## Parameters
 * The algorithm takes some parameters, which change the internal behaviour.
 * For the full list see the class GuessFromArchivedKeys
 * - "eigensolver":   Parameter submap which controls the behaviour of
 *                    the eigensolver. By default the solver will try to
 *                    find the smallest real eigenvalues
 *
 * \throws ExcObtainingScfGuessFailed if obtaining the guess failed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> guess_from_archived(
      const MolecularSystem& system, const IntegralOperator& fock_bb,
      const OverlapMatrix& S_bb, const krims::GenMap& params) {

  // Idea:   Do Sturmian calculation with minimal nmax as a dense matrix
  //         for a guess. Potentially even store the result *hard coded*
  //         in the program (That's how ORCA does it).
  //
  //         max(n_alpha+n_beta)    nmax
  //                          1     1
  //                     <=   5     2
  //                     <=  14     3
  //                     <=  30     4
  //                     <=  55     5
  //                     <=  91     6

  assert_throw(false,
               ExcObtainingScfGuessFailed("Method from_archived not yet implemented."));
  (void)system;
  (void)fock_bb;
  (void)S_bb;
  (void)params;

  return linalgwrap::EigensolutionTypeFor<true, IntegralOperator>{};
}

}  // namespace molsturm
