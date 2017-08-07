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
#include "RestrictionType.hh"
#include <gint/Integral.hh>
#include <krims/TypeUtils/UsingLibrary.hh>
#include <lazyten/BlockDiagonalMatrix.hh>

namespace molsturm {

// TODO  This class feels wrong and seems to indicate a bad and overly complicated design
//       in my opinion ....

/** Class representing the overlap matrix corresponding to a problem of the
 * provided RestrictionType */
template <typename StoredMatrix, RestrictionType restriction>
class OverlapMatrix {
  static_assert(restriction != RestrictionType::RestrictedOpen,
                "RestrictedOpen is not implemented yet.");
};

template <typename StoredMatrix>
class OverlapMatrix<StoredMatrix, RestrictionType::Unrestricted>
      : public lazyten::BlockDiagonalMatrix<gint::Integral<StoredMatrix>, 2> {
 public:
  typedef gint::Integral<StoredMatrix> int_term_type;
  OverlapMatrix(int_term_type S_bb)
        : lazyten::BlockDiagonalMatrix<int_term_type, 2>{{{S_bb, S_bb}}} {}
  const int_term_type& block_alpha() const { return this->diag_blocks()[0]; }
  const int_term_type& block_beta() const { return this->diag_blocks()[1]; }
};

template <typename StoredMatrix>
class OverlapMatrix<StoredMatrix, RestrictionType::RestrictedClosed>
      : public gint::Integral<StoredMatrix> {
 public:
  typedef gint::Integral<StoredMatrix> int_term_type;
  OverlapMatrix(int_term_type S_bb) : int_term_type(std::move(S_bb)) {}
  const int_term_type& block_alpha() const { return *this; }
  const int_term_type& block_beta() const { return block_alpha(); }
};

}  // namespace molsturm
