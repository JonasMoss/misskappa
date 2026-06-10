#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

using misskappa::IntMat;
using misskappa::na_code;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

// 10-subject 2-rater complete-data fixture (same as available_test.cpp).
IntMat ten_subject_2rater() {
  IntMat x(10, 2);
  x <<
    0, 0,
    0, 0,
    0, 1,
    1, 0,
    1, 1,
    1, 1,
    1, 1,
    0, 0,
    1, 1,
    0, 1;
  return x;
}

// 12-subject 3-rater 3-category fixture with two NA cells. Locked-in
// reference values come from the dev/legacy/misskappa build.
IntMat twelve_subject_3rater_3cat() {
  IntMat x(12, 3);
  // Mirrors the R fixture used to capture the reference values:
  //   xs <- cbind(rater1, rater2, rater3)
  //   xs[2, 3] <- NA_integer_   # row 2 (1-indexed), col 3
  //   xs[5, 1] <- NA_integer_   # row 5 (1-indexed), col 1
  x <<
    0, 0, 0,
    1, 1, na_code,
    2, 2, 2,
    0, 0, 1,
    na_code, 2, 1,
    2, 1, 2,
    0, 1, 0,
    1, 1, 2,
    2, 2, 2,
    0, 0, 0,
    1, 0, 1,
    2, 2, 1;
  return x;
}

}  // namespace

TEST_CASE("Complete data: available / IPW / Gwet agree exactly") {
  IntMat x = ten_subject_2rater();
  auto W = ms::loss::identity_weights(2);
  REQUIRE(W.has_value());

  auto r_av  = ms::estimate_available(x, *W);
  auto r_ipw = ms::estimate_ipw(x, *W);
  auto r_gw  = ms::estimate_gwet(x, *W);

  REQUIRE(r_av.has_value());
  REQUIRE(r_ipw.has_value());
  REQUIRE(r_gw.has_value());

  // With no missing entries, pi_j = 1 for all raters, so all three modes
  // collapse to the same estimator.
  CHECK((r_av->estimates - r_ipw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((r_av->estimates - r_gw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((r_av->vcov - r_ipw->vcov).cwiseAbs().maxCoeff() < tol);
  CHECK((r_av->vcov - r_gw->vcov).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("IPW with quadratic weights matches legacy on a 3-rater fixture") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  REQUIRE(W.has_value());

  auto r = ms::estimate_ipw(x, *W);
  REQUIRE(r.has_value());

  // Frozen against the dev/legacy/misskappa build (kappa_raw(method="ipw")).
  CHECK(std::abs(r->estimates(0) - 0.7022388059701492) < 1e-12);  // Conger
  CHECK(std::abs(r->estimates(1) - 0.7018342391304349) < 1e-12);  // Fleiss
  CHECK(std::abs(r->estimates(2) - 0.6977272727272728) < 1e-12);  // BP
}

TEST_CASE("Available with quadratic weights matches legacy on the same fixture") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);

  auto r = ms::estimate_available(x, *W);
  REQUIRE(r.has_value());

  CHECK(std::abs(r->estimates(0) - 0.6997960652591171) < 1e-12);
  CHECK(std::abs(r->estimates(1) - 0.6993437900128041) < 1e-12);
  CHECK(std::abs(r->estimates(2) - 0.6953125000000000) < 1e-12);
}

TEST_CASE("Gwet with quadratic weights matches legacy on the same fixture") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);

  auto r = ms::estimate_gwet(x, *W);
  REQUIRE(r.has_value());

  CHECK(std::abs(r->estimates(0) - 0.6998600746268656) < 1e-12);
  CHECK(std::abs(r->estimates(1) - 0.6994522758152174) < 1e-12);
  CHECK(std::abs(r->estimates(2) - 0.6953125000000000) < 1e-12);
}

TEST_CASE("IPW under uniform MCAR is close to available-case") {
  // 50-subject 3-rater dataset with uniform MCAR missingness across raters.
  // Under uniform MCAR the IPW and available-case estimators target the
  // same population quantity, so their estimates should agree to within
  // sampling noise at this sample size.
  IntMat x(50, 3);
  for (int i = 0; i < 50; ++i) {
    // Deterministic but mixed: ratings depend on i with a small "shift" per rater.
    x(i, 0) = i % 3;
    x(i, 1) = (i % 3 + (i / 3) % 2) % 3;
    x(i, 2) = (i % 3 + (i / 6) % 2) % 3;
  }
  // Uniform MCAR at 10% per cell, seeded for reproducibility.
  std::srand(42);
  for (int i = 0; i < 50; ++i) {
    for (int j = 0; j < 3; ++j) {
      if ((std::rand() % 10) == 0) x(i, j) = na_code;
    }
  }

  auto W = ms::loss::identity_weights(3);
  auto r_av = ms::estimate_available(x, *W);
  auto r_ipw = ms::estimate_ipw(x, *W);
  REQUIRE(r_av.has_value());
  REQUIRE(r_ipw.has_value());

  // Loose tolerance: same target population, different finite-sample paths.
  // The estimators are not numerically identical with non-uniform missingness
  // even if the design is MCAR (rater-specific pi_j are still random in finite
  // samples). 1e-2 is comfortable headroom at n=50, R=3.
  CHECK((r_av->estimates - r_ipw->estimates).cwiseAbs().maxCoeff() < 1e-2);
}

namespace {

void check_symmetric_psd(const RealMat& V) {
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-12);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-9);
}

}  // namespace

TEST_CASE("IPW variance is symmetric and PSD on quadratic-weight fixture") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_ipw(x, *W);
  REQUIRE(r.has_value());
  check_symmetric_psd(r->vcov);
}

TEST_CASE("Gwet variance is symmetric and PSD on quadratic-weight fixture") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_gwet(x, *W);
  REQUIRE(r.has_value());
  check_symmetric_psd(r->vcov);
}

TEST_CASE("IPW rejects a rater with zero observations") {
  IntMat x(5, 3);
  x << 0, 0, na_code,
       1, 1, na_code,
       0, 1, na_code,
       1, 0, na_code,
       0, 0, na_code;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_ipw(x, *W);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::not_identified);
}

TEST_CASE("Gwet rejects a rater with zero observations") {
  IntMat x(5, 3);
  x << 0, 0, na_code,
       1, 1, na_code,
       0, 1, na_code,
       1, 0, na_code,
       0, 0, na_code;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_gwet(x, *W);
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::not_identified);
}

TEST_CASE("Raw moment estimators reject incomplete rater co-observation graph") {
  IntMat x(6, 3);
  x << 0, 0, na_code,
       0, 1, na_code,
       1, 0, na_code,
       na_code, 0, 0,
       na_code, 0, 1,
       na_code, 1, 1;
  auto W = ms::loss::identity_weights(2);

  auto available = ms::estimate_available(x, *W);
  auto ipw = ms::estimate_ipw(x, *W);
  auto gwet = ms::estimate_gwet(x, *W);
  REQUIRE(!available.has_value());
  REQUIRE(!ipw.has_value());
  REQUIRE(!gwet.has_value());
  CHECK(available.error() == ms::Error::not_identified);
  CHECK(ipw.error() == ms::Error::not_identified);
  CHECK(gwet.error() == ms::Error::not_identified);
}
