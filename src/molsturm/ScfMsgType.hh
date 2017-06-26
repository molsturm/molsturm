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
#include <type_traits>

namespace molsturm {

/** How much and what should be printed during an SCF iteration */
enum class ScfMsgType {
  /** Print errors, but silence everything else */
  Silent = 0x0,
  /** Print some details summarising the iteration process */
  IterationProcess = 0x1,
  /** Print the final summary */
  FinalSummary = 0x2,
};

inline ScfMsgType operator&(ScfMsgType l, ScfMsgType r) {
  using type = typename std::underlying_type<ScfMsgType>::type;
  return static_cast<ScfMsgType>(static_cast<type>(l) & static_cast<type>(r));
}

inline ScfMsgType operator|(ScfMsgType l, ScfMsgType r) {
  using type = typename std::underlying_type<ScfMsgType>::type;
  return static_cast<ScfMsgType>(static_cast<type>(l) | static_cast<type>(r));
}

inline ScfMsgType& operator&=(ScfMsgType& l, ScfMsgType r) {
  l = l & r;
  return l;
}

inline ScfMsgType& operator|=(ScfMsgType& l, ScfMsgType r) {
  l = l | r;
  return l;
}

inline bool have_common_bit(ScfMsgType in, ScfMsgType what) {
  return static_cast<bool>(in & what);
}

}  // namespace molsturm
