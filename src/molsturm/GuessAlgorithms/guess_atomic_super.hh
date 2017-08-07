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

typedef GuessAlgorithmsKeysBase GuessAtomicSuperKeys;

/** Obtain an atomic superposition guess
 *
 * ## Parameters
 * The algorithm takes some parameters, which change the internal behaviour.
 * For the full list see the class GuessAtomicSuperKeys
 * - "eigensolver":   Parameter submap which controls the behaviour of
 *                    the eigensolver. By default the solver will try to
 *                    find the smallest real eigenvalues
 *
 * \throws ExcObtainingScfGuessFailed if obtaining the guess failed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
lazyten::EigensolutionTypeFor<true, IntegralOperator> guess_atomic_super(
      const MolecularSystem& system, const IntegralOperator& fock_bb,
      const OverlapMatrix& S_bb, const krims::GenMap& params) {

  assert_throw(false,
               ExcObtainingScfGuessFailed("Method atomic_super not yet implemented."));
  (void)system;
  (void)fock_bb;
  (void)S_bb;
  (void)params;

  return lazyten::EigensolutionTypeFor<true, IntegralOperator>{};
}

}  // namespace molsturm
