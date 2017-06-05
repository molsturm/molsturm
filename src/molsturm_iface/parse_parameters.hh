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
#include "ExcInvalidParameters.hh"
#include "Parameters.hh"
#include "molsturm/MolecularSystem.hh"
#include <gint/Structure.hh>
#include <krims/GenMap.hh>

namespace molsturm {
namespace iface {

/** Build the molecular system object from the parameters */
MolecularSystem build_molecular_system(const Parameters& params);

/** Examine the parameters and determine whether a restricted calculation
 *  shall be attempted or not */
bool parse_restricted(const Parameters& params, const MolecularSystem& system);

/** Build the integral lookup parameters from the input parameters and the molecular
 * system */
krims::GenMap build_int_params(const Parameters& params, const MolecularSystem& system);

/** Build the parameters for obtaining the scf guess */
krims::GenMap build_guess_params(const Parameters& params, const MolecularSystem& system);

/** Build the parameters for running the scf */
krims::GenMap build_scf_params(const Parameters& params, const MolecularSystem& system);

}  // namespace iface
}  // namespace molsturm
