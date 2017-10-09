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
#include "ScfParameters.hh"
#include "ScfResults.hh"
#include "ScfSolutionView.hh"

namespace molsturm {
namespace iface {

/** Run a SCF calculation */
ScfResults self_consistent_field(const ScfParameters& params,
                                 ScfSolutionView& solution_view);

}  // namespace iface
}  // namespace molsturm
