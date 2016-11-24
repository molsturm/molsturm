#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <krims/ExceptionSystem.hh>
#include <krims/NumComp.hh>

int main(int argc, char* const argv[]) {
  // Make sure that the program does not get aborted,
  // but all krims exceptions throw instead.
  krims::AssertDbgEffect::set_throw();

  // Throw in case a numerical comparison fails with very detailed
  // information
  krims::NumCompConstants::default_failure_action =
        krims::NumCompActionType::ThrowVerbose;

  // Run catch:
  int result = Catch::Session().run(argc, argv);
  return result;
}
