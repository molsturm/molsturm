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
#include <linalgwrap/ortho.hh>
#include <linalgwrap/random.hh>

namespace molsturm {

struct GuessRandomKeys {};

/** Obtain a random guess: Use random coefficient matrix
 *
 * Use a random coefficient matrix and make it S-orthogonal.
 *
 * ## Parameters
 * The algorithm takes currently no parameters.
 *
 * \throws ExcObtainingScfGuessFailed if obtaining the guess failed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> guess_random(
      const MolecularSystem& /*system*/, const IntegralOperator& fock_bb,
      const OverlapMatrix& S_bb, const krims::GenMap& /* params */) {
  using linalgwrap::random;
  using linalgwrap::ortho;

  typedef typename IntegralOperator::stored_matrix_type stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;

  // TODO Make alpha and beta block both truely random for unrestricted and not one
  // identical to the other
  const auto occa = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
  const auto occb = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA);
  const size_t n_vectors = std::max(occa.size(), occb.size());
  assert_implemented(occa == occb);

  // Get alpha-alpha block of the overlap matrix.
  // Note by construction this block is identical to the beta-beta block
  const auto& Sa_bb = S_bb.block_alpha();

  // Restricted open-shell is not yet implemented
  assert_internal(occa == occb || !fock_bb.restricted());

  linalgwrap::MultiVector<vector_type> guess;
  guess.reserve(n_vectors);
  for (size_t i = 0; i < n_vectors; ++i) {
    guess.push_back(random<vector_type>(Sa_bb.n_rows()));
  }

  linalgwrap::EigensolutionTypeFor<true, IntegralOperator> sol;
  sol.evectors() = ortho(guess, Sa_bb);
  sol.evalues() = std::vector<scalar_type>(n_vectors, 1);
  return fock_bb.restricted() ? sol : replicate_block(sol);
}

}  // namespace molsturm
