#include <gint/version.hh>
#include <gscf/version.hh>
#include <iostream>
#include <linalgwrap/version.hh>
#include <molsturm/version.hh>

namespace hard_coded_hf {}  // namespace hard_coded_hf

int main() {
  std::cout << "molsturm version:   " << molsturm::version::version_string()
            << std::endl
            << "gscf version:       " << gscf::version::version_string()
            << std::endl
            << "gint version:       " << gint::version::version_string()
            << std::endl
            << "linalgwrap version: " << linalgwrap::version::version_string()
            << std::endl;

  return 0;
}
