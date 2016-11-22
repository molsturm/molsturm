#pragma once
#include "molsturm/version_defs.hh"  // will be created by cmake
#include <string>

// The file molsturm/version_defs.hh will be created by cmake from
// the interal project version and will contain definitions of the
// macros VERSION_MINOR VERSION_MAJOR VERSION_PATCH

namespace molsturm {

struct version {
  static int constexpr major{detail::__version_var_major};
  static int constexpr minor{detail::__version_var_minor};
  static int constexpr patch{detail::__version_var_patch};

  // Return the version as a string
  static std::string version_string();
};

}  // namespace molsturm
