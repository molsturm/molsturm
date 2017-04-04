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
#include <gint/Structure.hh>

namespace molsturm {

template <typename Scalar = double>
Scalar energy_nuclear_repulsion(const gint::Structure& st) {
  // TODO
  // Some data structure for small vectors (e.g. points)
  // would make this easier to code.
  Scalar ret{0};
  for (auto itA = std::begin(st); itA != std::end(st); ++itA) {
    for (auto itB = itA + 1; itB != std::end(st); ++itB) {
      const std::array<Scalar, 3> diff{{itA->coords[0] - itB->coords[0],
                                        itA->coords[1] - itB->coords[1],
                                        itA->coords[2] - itB->coords[2]}};
      const Scalar rAB2 = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
      ret += itA->nuclear_charge * itB->nuclear_charge / std::sqrt(rAB2);
    }  // itB
  }    // itA
  return ret;
}

}  // namespace molsturm
