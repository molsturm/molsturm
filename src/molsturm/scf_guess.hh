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
#include "GuessAlgorithms.hh"

namespace molsturm {

struct ScfGuessKeys : public GuessAlgorithmsKeysBase {
  /** Guess method to use (Type: string) */
  static const std::string method;
};

/** Obtain a guess solution for an Scf Problem
 *
 * ## Control parameters and their default values
 * - method: The method to use for obtaining the guess.
 *   Avalable are:
 *      - "atomic_super":  Superposition of atomic densities
 *      - "extended_hueckel": Extended hueckel guess
 *      - "from_archived": Try to re-use an archived scf solution.
 *      - "hcore":    Diagonalise core hamiltonian
 *      - "loewdin":  Diagonalise S and use the eigenvectors
 *    (Default: hcore)
 * - eigensolver: Parameters for the eigensolver which is used.
 *
 * Some methods support even more parameters. See their
 * Keys object or their function documentation for details.
 */
template <typename IntegralOperator, typename OverlapMatrix>
linalgwrap::EigensolutionTypeFor<true, IntegralOperator> scf_guess(
      const MolecularSystem& system, const IntegralOperator& fock_bb,
      const OverlapMatrix& S_bb, const krims::GenMap& params = krims::GenMap{}) {
  const std::string method = params.at<std::string>(ScfGuessKeys::method, "hcore");

  if (method == std::string("atomic_super")) {
    return guess_atomic_super(system, fock_bb, S_bb, params);
  } else if (method == std::string("extended_hueckel")) {
    return guess_extended_hueckel(system, fock_bb, S_bb, params);
  } else if (method == std::string("from_archived")) {
    return guess_from_archived(system, fock_bb, S_bb, params);
  } else if (method == std::string("hcore")) {
    return guess_hcore(system, fock_bb, S_bb, params);
  } else if (method == std::string("loewdin")) {
    return guess_loewdin(system, fock_bb, S_bb, params);
  } else if (method == std::string("random")) {
    return guess_random(system, fock_bb, S_bb, params);
  } else {
    assert_throw(false, ExcInvalidScfGuessParametersEncountered(
                              "The method '" + method +
                              "' is not known. Did you spell it wrong?"));
    return linalgwrap::EigensolutionTypeFor<true, IntegralOperator>{};
  }
}

}  // namespace molsturm
