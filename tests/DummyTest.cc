#include <catch.hpp>
#include <rapidcheck.h>

namespace molsturm {
namespace tests {

TEST_CASE("Dummy test", "[dummy]") {
  auto test = [](int x) { RC_ASSERT(x + x == 2 * x); };
  CHECK(rc::check("Run dummy", test));
}

}  // namespace tests
}  // namespace molsturm
