#include "molsturm/version.hh"
#include <sstream>

namespace molsturm {

std::string version::version_string() {
  std::stringstream ss;
  ss << major << "." << minor << "." << patch;
  return ss.str();
}

}  // namespace molsturm
