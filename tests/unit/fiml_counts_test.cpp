#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

using misskappa::EmOptions;
using misskappa::IntMat;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

IntMat fleiss1971() {
  IntMat x(30, 5);
  x <<
    0, 0, 0, 6, 0,
    0, 3, 0, 0, 3,
    0, 1, 4, 0, 1,
    0, 0, 0, 0, 6,
    0, 3, 0, 3, 0,
    2, 0, 4, 0, 0,
    0, 0, 4, 0, 2,
    2, 0, 3, 1, 0,
    2, 0, 0, 4, 0,
    0, 0, 0, 0, 6,
    1, 0, 0, 5, 0,
    1, 1, 0, 4, 0,
    0, 3, 3, 0, 0,
    1, 0, 0, 5, 0,
    0, 2, 0, 3, 1,
    0, 0, 5, 0, 1,
    3, 0, 0, 1, 2,
    5, 1, 0, 0, 0,
    0, 2, 0, 4, 0,
    1, 0, 2, 0, 3,
    0, 0, 0, 0, 6,
    0, 1, 0, 5, 0,
    0, 2, 0, 1, 3,
    2, 0, 0, 4, 0,
    1, 0, 0, 4, 1,
    0, 5, 0, 1, 0,
    4, 0, 0, 0, 2,
    0, 2, 0, 4, 0,
    1, 0, 5, 0, 0,
    0, 0, 0, 0, 6;
  return x;
}

// 8 subjects, 3 categories, r_total = 4 with varying r_i so partial-data
// EM weighting actually matters.
IntMat partial_counts_fixture() {
  IntMat y(8, 3);
  y <<
    3, 1, 0,
    0, 2, 1,
    2, 0, 1,
    1, 2, 0,
    0, 0, 3,
    2, 1, 0,
    3, 0, 0,
    0, 1, 2;
  return y;
}

}  // namespace

TEST_CASE("FIML-counts on Fleiss 1971 (complete) matches available-case bit-for-bit") {
  IntMat x = fleiss1971();
  auto W = ms::loss::identity_weights(5);
  REQUIRE(W.has_value());

  // r_total = 6 for Fleiss 1971; every row sums to 6 so no partial counts.
  auto r = ms::estimate_fiml_counts(x, *W, /*r_total=*/6, EmOptions{});
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 2);

  // Frozen against dev/legacy/misskappa kappa_counts(method="ml").
  CHECK(std::abs(r->estimates(0) - 0.4302445200601408) < tol);  // Fleiss
  CHECK(std::abs(r->estimates(1) - 0.4444444444444444) < tol);  // BP

  // With complete counts FIML's vcov coincides with the moment-based one.
  CHECK(std::abs(r->vcov(0, 0) - 0.002839607123962024) < 1e-9);
  CHECK(std::abs(r->vcov(1, 1) - 0.002937242798353910) < 1e-9);
}

TEST_CASE("FIML-counts on Fleiss 1971 (complete) matches available with quadratic weights") {
  IntMat x = fleiss1971();
  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto W = ms::loss::quadratic_weights(5, v);

  auto r = ms::estimate_fiml_counts(x, *W, /*r_total=*/6, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.2840722495894914) < 1e-8);  // Fleiss
  CHECK(std::abs(r->estimates(1) - 0.3338888888888891) < 1e-8);  // BP
}

TEST_CASE("FIML-counts on a partial-data fixture differs from available-case") {
  IntMat y = partial_counts_fixture();
  auto W = ms::loss::identity_weights(3);

  auto r_fiml = ms::estimate_fiml_counts(y, *W, /*r_total=*/4, EmOptions{});
  REQUIRE(r_fiml.has_value());

  // Frozen against dev/legacy/misskappa kappa_counts(method="ml").
  CHECK(std::abs(r_fiml->estimates(0) - 0.2725053215314407) < 1e-7);
  CHECK(std::abs(r_fiml->estimates(1) - 0.2812687432376064) < 1e-7);

  // Variance matches legacy to a looser tolerance (Louis IF sensitive to
  // pruning + tol).
  CHECK(std::abs(r_fiml->vcov(0, 0) - 0.02307528179093083) < 1e-5);
  CHECK(std::abs(r_fiml->vcov(1, 1) - 0.02234091734312603) < 1e-5);
}

TEST_CASE("FIML-counts: variance is symmetric and PSD") {
  IntMat y = partial_counts_fixture();
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_fiml_counts(y, *W, /*r_total=*/4, EmOptions{});
  REQUIRE(r.has_value());
  const RealMat& V = r->vcov;
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-10);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-8);
}

TEST_CASE("FIML-counts: info_rcond affects Louis variance, not estimates") {
  IntMat y = partial_counts_fixture();
  auto W = ms::loss::identity_weights(3);

  EmOptions keep{};
  keep.info_rcond = 0.0;
  auto r_keep = ms::estimate_fiml_counts(y, *W, /*r_total=*/4, keep);
  REQUIRE(r_keep.has_value());

  EmOptions drop{};
  drop.info_rcond = 1.0;
  auto r_drop = ms::estimate_fiml_counts(y, *W, /*r_total=*/4, drop);
  REQUIRE(r_drop.has_value());

  CHECK((r_keep->estimates - r_drop->estimates).cwiseAbs().maxCoeff() < 1e-12);
  CHECK(r_drop->vcov.norm() < r_keep->vcov.norm());
}

TEST_CASE("FIML-counts: perfect agreement -> kappa = 1") {
  IntMat x(4, 3);
  x << 4, 0, 0,
       0, 4, 0,
       0, 0, 4,
       4, 0, 0;
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_fiml_counts(x, *W, /*r_total=*/4, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < 1e-6);  // Fleiss
  CHECK(std::abs(r->estimates(1) - 1.0) < 1e-6);  // BP
}

TEST_CASE("FIML-counts: dimension mismatch -> dimension_mismatch") {
  IntMat x(4, 3);
  x << 4, 0, 0,
       0, 4, 0,
       0, 0, 4,
       2, 2, 0;
  auto W = ms::loss::identity_weights(4);  // mismatched
  auto r = ms::estimate_fiml_counts(x, *W, /*r_total=*/4, EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::dimension_mismatch);
}

TEST_CASE("FIML-counts: row sum exceeding r_total -> invalid_argument") {
  IntMat x(3, 3);
  x << 4, 0, 0,
       0, 4, 0,
       5, 0, 0;  // row sums to 5 > r_total = 4
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_fiml_counts(x, *W, /*r_total=*/4, EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("FIML-counts: r_total < 2 -> invalid_argument") {
  IntMat x(2, 3);
  x << 1, 0, 0,
       0, 1, 0;
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_fiml_counts(x, *W, /*r_total=*/1, EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}
