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
#include <linalgwrap/SmallVector.hh>

namespace molsturm {

struct GuessExternalKeys : public GuessAlgorithmsKeysBase {
  /** The external eigensolution guess.
   *  Type: linalgwrap::Eigensolution<double, SmallVector<double>>.
   */
  const static std::string eigensolution;
};

/** Obtain a guess by using an externally provided guess solution.
 *
 * The solution is not altered. The exception is if the guess only holds an
 * alpha block, but we are supposed to provide a guess for an unrestricted
 * calculation. In this case the solution is duplicated onto the beta block.
 *
 * ## Parameters
 * The algorithm takes some parameters, which change the internal behaviour.
 * For the full list see the class GuessAlgorithmsKeysBase.
 *
 * - "eigensolver":   Unused.
 * - "eigensolution": The external eigensolution, which is to be used.
 *                    Mandatory parameter.
 *
 * \throws ExcObtainingScfGuessFailed if obtaining the guess failed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> guess_external(
      const MolecularSystem& /*system*/, const IntegralOperator& fock_bb,
      const OverlapMatrix& /*S_bb*/, const krims::GenMap& params) {
  typedef typename IntegralOperator::stored_matrix_type stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef typename stored_matrix_type::scalar_type scalar_type;
  static_assert(std::is_same<scalar_type, double>::value,
                "Scalar type needs to be double for guess_external at the moment.");

  assert_throw(params.exists(GuessExternalKeys::eigensolution),
               ExcInvalidScfGuessParametersEncountered(
                     "For guess_external the parameter eigensolution is mandatory."));

  // Check eigensolution is as expected
  auto esolution =
        params.at<linalgwrap::Eigensolution<double, linalgwrap::SmallVector<double>>>(
              GuessExternalKeys::eigensolution);

  auto& evectors = esolution.evectors();
  auto& evalues = esolution.evalues();
  assert_throw(evectors.n_vectors() == evalues.size(),
               ExcInvalidScfGuessParametersEncountered(
                     "Number of eigenvectors and number of eigenvalues in provided "
                     "eigensolution does not agree."));

  const auto occa = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
  const auto occb = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA);
  const size_t n_vectors = std::max(occa.length(), occb.length());

  // Restricted open-shell is not yet implemented
  assert_internal(occa == occb || !fock_bb.restricted());
  if (fock_bb.restricted()) {
    assert_throw(
          evectors.n_vectors() == n_vectors,
          ExcInvalidScfGuessParametersEncountered("For restricted the number of "
                                                  "eigenvectors needs to be "
                                                  "exactly max(n_alpha,n_beta) == " +
                                                  std::to_string(n_vectors)));
  } else {
    assert_throw(
          evectors.n_vectors() == n_vectors || evectors.n_vectors() == 2 * n_vectors,
          ExcInvalidScfGuessParametersEncountered(
                "For unrestricted the number of eigenvectors either needs to be "
                "max(n_alpha,n_beta) == " +
                std::to_string(n_vectors) +
                " or twice this value. In the former case the "
                "guess will be duplicated into the beta block."));
  }

  // Copy into the datastructure we use and replicate the block if needed
  linalgwrap::EigensolutionTypeFor<true, IntegralOperator> sol;
  sol.evalues() = evalues;
  sol.evectors().reserve(evalues.size());
  for (auto& vec : evectors) {
    vector_type nvec(vec.size());
    std::copy(vec.begin(), vec.end(), nvec.begin());
    sol.evectors().push_back(std::move(nvec));
  }

  if (fock_bb.restricted() || sol.evectors().n_vectors() == 2 * n_vectors) {
    return sol;
  } else {
    return replicate_block(sol);
  }
}

}  // namespace molsturm
