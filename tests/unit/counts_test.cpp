#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

using misskappa::IntMat;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

// 8-subject 3-category count fixture; row sums are 4 for every subject
// (so r_i = 4, every subject rated by 4 raters).
IntMat small_counts() {
  IntMat x(8, 3);
  x << 4, 0, 0,
       0, 4, 0,
       0, 0, 4,
       3, 1, 0,
       1, 3, 0,
       0, 3, 1,
       2, 1, 1,
       1, 2, 1;
  return x;
}

}  // namespace

TEST_CASE("Counts available-case matches legacy on Fleiss-shaped data") {
  // Mirror of the Fleiss 1971 5-category counts. Hand-typed from
  // dev/legacy/misskappa/data/dat.fleiss1971.rda; subjects are rows,
  // categories columns 1..5, row sums = 6 raters per subject.
  IntMat x(30, 5);
  // Verbatim from dev/legacy/misskappa/data/dat.fleiss1971; 30 subjects,
  // 6 raters, 5 diagnostic categories.
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

  auto W = ms::loss::identity_weights(5);
  REQUIRE(W.has_value());

  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 2);

  // Frozen against dev/legacy/misskappa kappa_counts(method="available",
  // weight="identity"). Classic Fleiss-1971 value: Fleiss ~= 0.4302.
  CHECK(std::abs(r->estimates(0) - 0.4302445200601408) < tol);  // Fleiss
  CHECK(std::abs(r->estimates(1) - 0.4444444444444444) < tol);  // BP

  CHECK(std::abs(r->vcov(0, 0) - 0.002839607123962025) < 1e-12);
  CHECK(std::abs(r->vcov(1, 1) - 0.002937242798353909) < 1e-12);
  CHECK(std::abs(r->vcov(0, 1) - 0.002826379113973653) < 1e-12);
}

TEST_CASE("Counts available-case: quadratic weights match legacy") {
  IntMat x(30, 5);
  // Verbatim from dev/legacy/misskappa/data/dat.fleiss1971; 30 subjects,
  // 6 raters, 5 diagnostic categories.
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

  RealVec v(5);
  v << 1.0, 2.0, 3.0, 4.0, 5.0;
  auto W = ms::loss::quadratic_weights(5, v);

  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.2840722495894910) < tol);
  CHECK(std::abs(r->estimates(1) - 0.3338888888888889) < tol);
}

TEST_CASE("Counts available-case: perfect agreement -> kappa = 1") {
  // Every subject's counts are concentrated on one category (perfect rater
  // agreement). All rater pairs agree, so kappa = 1.
  IntMat x(4, 3);
  x << 4, 0, 0,
       0, 4, 0,
       0, 0, 4,
       4, 0, 0;
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);  // Fleiss
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);  // BP
}

TEST_CASE("Counts available-case: variance is symmetric and PSD") {
  IntMat x = small_counts();
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(r.has_value());
  const RealMat& V = r->vcov;
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-12);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-10);

  REQUIRE(r->psi.rows() == x.rows());
  REQUIRE(r->psi.cols() == 2);
  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - V).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("Counts: dimension mismatch -> dimension_mismatch") {
  IntMat x = small_counts();
  // W is 4x4 but x has 3 columns.
  auto W = ms::loss::identity_weights(4);
  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::dimension_mismatch);
}

TEST_CASE("Counts: negative entry -> invalid_argument") {
  IntMat x(3, 3);
  x << 2, 1, 1,
       1, -1, 2,
       0, 3, 1;
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("Counts: varying row sums work (partial counts)") {
  IntMat x(3, 3);
  x << 4, 0, 0,
       2, 1, 0,  // r_i = 3
       1, 1, 1;  // r_i = 3
  auto W = ms::loss::identity_weights(3);
  auto r = ms::estimate_available_counts(x, *W);
  REQUIRE(r.has_value());
  for (Eigen::Index i = 0; i < r->estimates.size(); ++i) {
    CHECK(std::isfinite(r->estimates(i)));
  }
}
