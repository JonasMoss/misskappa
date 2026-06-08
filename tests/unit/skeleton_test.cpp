// Skeleton test: validates that the no-exceptions toolchain compiles a
// doctest binary that links the public surface. Replaced by real loss /
// estimator tests as those modules come online.

#include "doctest.h"

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

TEST_CASE("Result<int> carries values") {
  misskappa::Result<int> ok = 42;
  REQUIRE(ok.has_value());
  CHECK(*ok == 42);
}

TEST_CASE("Result<int> carries errors") {
  misskappa::Result<int> bad = misskappa::unexpected(misskappa::Error::invalid_argument);
  REQUIRE(!bad.has_value());
  CHECK(bad.error() == misskappa::Error::invalid_argument);
}

TEST_CASE("na_code is -1") {
  CHECK(misskappa::na_code == -1);
}
