#pragma once
#include "molsturm/config.hh"
#include <string>

namespace molsturm {

struct Version {
  static int const major{detail::version_major};
  static int const minor{detail::version_minor};
  static int const patch{detail::version_patch};

  // Return the version as a string
  static std::string as_string();
};

}  // namespace molsturm
