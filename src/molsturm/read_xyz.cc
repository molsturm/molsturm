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

#include "read_xyz.hh"
#include "gint/Element.hh"

namespace molsturm {

gint::Structure read_xyz(std::istream& f, double angstrom_to_bohr) {
  static_assert(std::is_same<gint::real_type, double>::value, "Type mismatch error");
  assert_throw(f, krims::ExcIO());

  size_t n_atoms;
  f >> n_atoms;
  assert_throw(f,
               ExcInvalidXyz("First line cannot be interpreted as the number of atoms"));

  // The second line is a comment => read and ignore
  std::string comment;
  std::getline(f, comment);

  // Now come the atoms:
  gint::Structure molec;
  molec.reserve(n_atoms);

  for (size_t i = 0; i < n_atoms; ++i) {
    std::string symbol;
    double x, y, z;
    f >> symbol >> x >> y >> z;
    assert_throw(f, ExcInvalidXyz("Error at line " + std::to_string(2 + i) + " of atom " +
                                  symbol + "."));

    if (gint::is_element_symbol(symbol)) {
      // Convert units and place into molecule
      molec.push_back(gint::Atom(symbol, {{x * angstrom_to_bohr, y * angstrom_to_bohr,
                                           z * angstrom_to_bohr}}));
    } else {
      assert_throw(false, ExcInvalidXyz("Error at line: Unknown element symbol \"" +
                                        symbol + "\"."));
    }
  }

  assert_throw(
        n_atoms == molec.size(),
        ExcInvalidXyz(
              "Numbor of atoms read and number of atoms to be read does not agree."));
  return molec;
}
}  // namespace molsturm
