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
