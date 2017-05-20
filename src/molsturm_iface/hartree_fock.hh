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
#include "HfResults.hh"
#include "Parameters.hh"
#include <krims/ExceptionSystem.hh>

namespace molsturm {
namespace iface {

#ifndef SWIG
DefException2(ExcTooSmallBasis, size_t, size_t,
              << "A basis of size " << arg1
              << " is too small to incorporate max(alpha,beta) = " << arg2
              << " electrons. Choose a larger basis.");

DefException1(ExcInvalidParameters, std::string,
              << "Invalid parameters passed to molsturm: " << arg1);
#endif  // SWIG

HfResults hartree_fock(const Parameters& p);

}  // namespace iface
}  // namespace molsturm