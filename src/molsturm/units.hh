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

namespace molsturm {

namespace units {
// Taken from http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
// DOI: 10.5281/zenodo.22826
static constexpr double bohr_radius      = 0.52917721067e-10;
static constexpr double bohr_to_angstrom = 1e10 * bohr_radius;
static constexpr double angstrom_to_bohr = 1 / bohr_to_angstrom;

// http://physics.nist.gov/cgi-bin/cuu/Value?threv
// DOI: 10.5281/zenodo.22826
static constexpr double hartree_energy = 27.21138602;  // eV
static constexpr double hartree_to_eV  = hartree_energy;
static constexpr double eV_to_hartree  = 1 / hartree_to_eV;

}  // namespace units

}  // namespace molsturm
