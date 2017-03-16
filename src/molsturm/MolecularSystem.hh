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
#include <gint/chemistry/Structure.hh>
#include <linalgwrap/Constants.hh>

namespace molsturm {

DefException2(ExcInvalidMultiplicity, size_t, std::string,
              << "The multiplicity \"" << arg1
              << " is invalid for this molecular system: " << arg2);

DefException1(ExcNonIntegerElectronCount, double, << "The determined electron count "
                                                  << arg1 << " is not an integer value.");

// TODO Allow to construct from GenMap!
/** Class which describes the physical system to model
 *  in this calculation */
struct MolecularSystem {
  /** The structure, i.e. the arrangement of the nuclei */
  gint::Structure structure;

  /** The number of alpha electrons */
  size_t n_alpha;

  /** The number of beta electrons */
  size_t n_beta;

  /** The total charge of the system */
  double charge;

  /** Construct an empty system, i.e. one with no electrons, charge or atoms */
  MolecularSystem() : structure{}, n_alpha(0), n_beta(0), charge(0) {}

  /** Construct a system from a structure and the number of alpha and beta electrons
   *
   * \param n_electrons  The number of electrons, first the number of alphas than
   *                     the number of betas. If number of betas is missing,
   *                     the number of alphas is used.
   **/
  //@{
  MolecularSystem(gint::Structure structure_, std::array<size_t, 2> n_electrons)
        : structure(std::move(structure_)),
          n_alpha(n_electrons[0]),
          n_beta(n_electrons[1]),
          charge(structure.total_charge() - n_electrons[0] - n_electrons[1]) {}
  //@}

  /** Construct a system from a structure, the total charge and the multiplicity
   *
   * If the multiplicity is missing 1 is assumed for systems with an even number
   * of electrons and 2 is assumed for a system with an odd number of electrons.
   *
   * If the charge is missing 0 is assumed.
   */
  explicit MolecularSystem(gint::Structure structure_, double charge_ = 0.,
                           size_t multiplicity_ = linalgwrap::Constants<size_t>::invalid);
};

}  // namespace molsturm
