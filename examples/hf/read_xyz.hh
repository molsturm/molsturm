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
#include <gint/chemistry/Molecule.hh>
#include <iostream>
#include <krims/ExceptionSystem.hh>

namespace hf {

// TODO put this into a namespace in molsturm?
namespace units {

// Taken from http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
static constexpr double bohr_radius = 0.52917721067e-10;
static constexpr double bohr_to_angstrom = 1e10 * bohr_radius;
static constexpr double angstrom_to_bohr = 1 / bohr_to_angstrom;
}

DefException1(ExcInvalidXyz, std::string,
              << "The input stream does not contain a valid xyz file: " << arg1);

/** Read a xyz file from a file stream
 *
 * \param bohr_to_angstrom   The conversion factor to use for converting the xyz file
 *                           coordinates (by convention: Angström) to bohr.
 *                           Set this to 1 if you use atomic units (bohr) in your input
 *                           file.
 *
 * */
gint::Molecule read_xyz(std::istream& f,
                        double angstrom_to_bohr = units::angstrom_to_bohr);

}  // namespace hf
