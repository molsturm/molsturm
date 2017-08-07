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

struct GuessHcoreKeys : public GuessAlgorithmsKeysBase {
  /** Prefactor to use in front of the nuclear attraction operator
   *  (type: scalar_type) */
  static const std::string nuclear_attraction_factor;
};

/** Obtain a core Hamiltonian guess
 *
 * ## Parameters
 * The algorithm takes some parameters, which change the internal behaviour.
 * For the full list see the class GuessHcoreKeys
 * - "eigensolver":   Parameter submap which controls the behaviour of
 *                    the eigensolver. By default the solver will try to
 *                    find the smallest real eigenvalues
 * - "nuclear_attraction_factor": The prefactor to use in front of the
 *                    nuclear attraction operator in the hcore guess.
 *
 * \throws ExcObtainingScfGuessFailed if obtaining the guess failed.
 */
template <typename IntegralOperator, typename OverlapMatrix>
lazyten::EigensolutionTypeFor<true, IntegralOperator> guess_hcore(
      const MolecularSystem& /*system*/, const IntegralOperator& fock_bb,
      const OverlapMatrix& S_bb, const krims::GenMap& params) {
  using namespace lazyten;
  typedef typename OverlapMatrix::scalar_type scalar_type;

  // TODO We had some trouble when choosing the factor different from 1.0,
  //      where the SCF all of a sudden started to converge to a wrong
  //      ground state, so this is changed back to 1.0 for now.
  const scalar_type nuc_attr_fac =
        params.at<scalar_type>(GuessHcoreKeys::nuclear_attraction_factor, 1.0);

  krims::GenMap eigensolver_params = params.submap(GuessHcoreKeys::eigensolver_params);
  eigensolver_params.insert_default(EigensystemSolverKeys::which, "SR");

  // Setup Hcore (i.e. only take one electron terms from Fock)
  // Scale the nuclear attraction term according to the
  LazyMatrixSum<typename IntegralOperator::stored_matrix_type> hcore;
  for (const auto& id_term : fock_bb.terms_1e()) {
    if (id_term.first.integral_type() == gint::IntegralType::nuclear_attraction) {
      hcore += nuc_attr_fac * id_term.second;
    } else {
      hcore += id_term.second;
    }
  }

  const auto occa        = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_ALPHA);
  const auto occb        = fock_bb.indices_orbspace(gscf::OrbitalSpace::OCC_BETA);
  const size_t n_vectors = std::max(occa.length(), occb.length());

  // Get alpha-alpha block of the overlap matrix.
  // Note by construction this block is identical to the beta-beta block
  const auto& Sa_bb = S_bb.block_alpha();

  // Restricted open-shell is not yet implemented
  assert_internal(occa == occb || !fock_bb.restricted());
  try {
    auto sol = eigensystem_hermitian(hcore, Sa_bb, n_vectors, eigensolver_params);
    return fock_bb.restricted() ? sol : replicate_block(sol);
  } catch (const SolverException& e) {
    rescue::failed_eigenproblem(
          Eigenproblem<true, decltype(hcore), typename OverlapMatrix::int_term_type>(
                hcore, Sa_bb),
          params);
    assert_throw(false, ExcObtainingScfGuessFailed(
                              "Eigensolver for Hcore failed with message " + e.extra()));
    return EigensolutionTypeFor<true, IntegralOperator>{};
  }
}

}  // namespace molsturm
