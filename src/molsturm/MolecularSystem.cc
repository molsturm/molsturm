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

#include "MolecularSystem.hh"

namespace molsturm {

MolecularSystem::MolecularSystem(gint::Structure structure_, double charge_,
                                 size_t multiplicity_)
      : structure(std::move(structure_)),
        n_alpha(0),
        n_beta(0),
        charge(std::move(charge_)) {
  // Compute number of electrons to distribute to alpha and beta spin:
  const double n_elec = structure.total_charge() - charge;
  const size_t n_elec_count = static_cast<size_t>(n_elec);
  assert_throw(n_elec - n_elec_count < 1e-12, ExcNonIntegerElectronCount(n_elec));

  // Twice the total spin (i.e. the multiplicity reduced by one)
  const size_t spin_twice = [&multiplicity_, &n_elec_count]() {
    if (multiplicity_ == linalgwrap::Constants<size_t>::invalid) {
      // Set to 0 if even number of electrons and 1 if odd
      // equals multiplicities 1 and 2, respectively
      return n_elec_count % 2 == 0 ? size_t(0) : size_t(1);
    } else {
      return multiplicity_ - 1;
    }
  }();

  assert_throw(
        multiplicity_ == linalgwrap::Constants<size_t>::invalid || multiplicity_ > 0,
        ExcInvalidMultiplicity(multiplicity_,
                               "The multiplicity needs to be a positive number."));
  assert_throw(
        spin_twice <= n_elec_count,
        ExcInvalidMultiplicity(
              multiplicity_,
              "The multiplicity should be larger or equal to the electron count (" +
                    std::to_string(n_elec_count) + ") plus 1"));

  const bool even_elec = n_elec_count % 2 == 0;
  const bool even_2spin = spin_twice % 2 == 0;

  // Odd spin_twice values imply a half-integer total spin
  // In other words we need an odd number of electrons in this case
  //
  // Even spin_twice values imply an integer total spin
  // i.e. an even number of electrons
  assert_throw(even_2spin == even_elec,
               ExcInvalidMultiplicity(
                     multiplicity_,
                     "Only a system with an even number of electrons can have a odd "
                     "multiplicity and only a system with an odd number of electrons "
                     "can have an even multiplicity. This system has " +
                           std::to_string(n_elec_count) + " electrons."));

  // Compute number of alpha and beta electrons
  n_alpha = (n_elec_count - spin_twice) / 2 + spin_twice;
  n_beta = n_elec_count - n_alpha;

  // Check that this makes sense
  if (multiplicity_ == linalgwrap::Constants<size_t>::invalid) {
    assert_internal((n_alpha - n_beta) == spin_twice);
  } else {
    assert_internal((n_alpha - n_beta) + 1 == multiplicity_);
  }
}

}  // namespace molsturm
