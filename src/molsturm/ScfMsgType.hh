#pragma once

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

ScfMsgType operator&(ScfMsgType l, ScfMsgType r) {
  return static_cast<ScfMsgType>(static_cast<int>(l) & static_cast<int>(r));
}

ScfMsgType operator|(ScfMsgType l, ScfMsgType r) {
  return static_cast<ScfMsgType>(static_cast<int>(l) | static_cast<int>(r));
}

ScfMsgType& operator&=(ScfMsgType& l, ScfMsgType r) {
  l = l & r;
  return l;
}

ScfMsgType& operator|=(ScfMsgType& l, ScfMsgType r) {
  l = l | r;
  return l;
}

bool common_bit(ScfMsgType in, ScfMsgType what) { return static_cast<bool>(in & what); }

}  // namespace molsturm
