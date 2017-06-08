#include "molsturm/Version.hh"
#include <sstream>

namespace molsturm {

std::string Version::as_string() {
  std::stringstream ss;
  ss << major << "." << minor << "." << patch;
  return ss.str();
}

}  // namespace molsturm
