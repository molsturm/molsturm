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
#include "molsturm/units.hh"
#include <gint/Structure.hh>
#include <iostream>
#include <krims/ExceptionSystem.hh>

namespace molsturm {

DefException1(
      ExcInvalidXyz, std::string,
      << "The input stream provided to read_xyz does not contain a valid xyz file: "
      << arg1);

/** Read a xyz file from a file stream
 *
 * \param bohr_to_angstrom   The conversion factor to use for converting the xyz file
 *                           coordinates (by convention: Ångström) to bohr.
 *                           Set this to 1 if you use atomic units (bohr) in your
 *                           input file.
 * */
gint::Structure read_xyz(std::istream& f,
                         double angstrom_to_bohr = units::angstrom_to_bohr);

}  // namespace molsturm
