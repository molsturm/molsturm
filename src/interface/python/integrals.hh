//
// Copyright (C) 2018 by the molsturm authors
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

#ifdef SWIG
// clang-format off
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* out_bb, int n_bas1, int n_bas_2)};
// clang-format on
#endif  // SWIG
namespace molsturm {
namespace iface {

/** Get the overlap matrix in basis functions */
void overlap_bb(const ScfParameters& params, double* out_bb, int n_bas1, int n_bas2);

/** Get the kinetic energy matrix in basis functions */
void kinetic_bb(const ScfParameters& params, double* out_bb, int n_bas1, int n_bas2);

/** Get the nuclear attraction matrix in basis functions */
void nuclear_attraction_bb(const ScfParameters& params, double* out_bb, int n_bas1,
                           int n_bas2);

}  // namespace iface
}  // namespace molsturm
