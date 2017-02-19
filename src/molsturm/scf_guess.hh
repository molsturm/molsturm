//
// Copyright (C) 2016 by the molsturm authors
//
// This file is part of molsturm.
//
// molsturm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// molsturm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with molsturm. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "GuessLibrary.hh"

namespace molsturm {

/** Thrown if solver parameters are wrong or inconsistent */
DefException1(ExcInvalidScfGuessParametersEncountered, std::string,
              << "The Scf guess could not be obtained, because the set of "
                 "parameters passed is not valid. Details: "
              << arg1);

struct ScfGuessKeys {
  /** Guess method to use (Type: string) */
  static const std::string method;
};

/** Obtain a guess solution for an Scf Problem
 *
 * ## Control parameters and their default values
 * - method: The method to use for obtaining the guess.
 *   Avalable are:
 *      - "loewdin":  Diagonalise S and use the eigenvectors
 *      - "hcore":    Diagonalise core hamiltonian
 *    (Default: hcore)
 */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> scf_guess(
      const IntegralOperator& fock_bb, const OverlapMatrix& s_bb,
      const krims::GenMap& params) {
  const std::string method = params.at<std::string>(ScfGuessKeys::method, "hcore");

  if (method == std::string("loewdin")) {
    return loewdin_guess(fock_bb, s_bb);
  } else if (method == std::string("hcore")) {
    return hcore_guess(fock_bb, s_bb);
  } else {
    assert_throw(false, ExcInvalidScfGuessParametersEncountered(
                              "The method '" + method +
                              "' is not known. Did you spell it wrong?"));
  }
}

}  // namespace molsturm
