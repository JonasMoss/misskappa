#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>

using misskappa::EmOptions;
using misskappa::IntMat;
using misskappa::na_code;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

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

IntMat twelve_subject_3rater_3cat() {
  IntMat x(12, 3);
  // Mirrors the same R fixture as ipw_gwet_test.cpp:
  //   xs[2, 3] <- NA;  xs[5, 1] <- NA  (1-indexed positions).
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

IntMat frechet_fixture_int() {
  IntMat x(4, 5);
  x << 0, 0, 1, 0, 0,
       0, 1, 2, 1, 1,
       1, 0, 0, 0, 0,
       1, 2, 3, 3, 4;
  return x;
}

}  // namespace

TEST_CASE("estimate_fiml: complete data, identity weights matches Cohen") {
  IntMat x = ten_subject_2rater();
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());

  // Frozen against legacy: complete-data EM converges to the empirical
  // pattern frequencies, so FIML reproduces the moment-based kappa values.
  CHECK(std::abs(r->estimates(0) - 0.4) < 1e-9);                // Conger
  CHECK(std::abs(r->estimates(1) - 13.0 / 33.0) < 1e-9);        // Fleiss
  CHECK(std::abs(r->estimates(2) - 0.4) < 1e-9);                // BP
}

TEST_CASE("estimate_fiml: 3-rater MAR fixture matches legacy") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());

  // Frozen against the dev/legacy/misskappa kappa_raw(method="ml") build.
  CHECK(std::abs(r->estimates(0) - 0.6778523464560615) < 1e-6);  // Conger
  CHECK(std::abs(r->estimates(1) - 0.6767676743661067) < 1e-6);  // Fleiss
  CHECK(std::abs(r->estimates(2) - 0.6666666647950232) < 1e-6);  // BP

  // Diagonal variance entries should match to a looser tolerance (EM SE is
  // sensitive to the pruning threshold and tol settings).
  CHECK(std::abs(r->vcov(0, 0) - 0.02971461597182125) < 1e-5);
  CHECK(std::abs(r->vcov(1, 1) - 0.02976059214783952) < 1e-5);
  CHECK(std::abs(r->vcov(2, 2) - 0.02387152699792634) < 1e-5);
}

TEST_CASE("estimate_fiml: variance is symmetric and PSD") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());
  const RealMat& V = r->vcov;
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-10);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-8);
}

TEST_CASE("estimate_fiml: influence functions reconstruct Louis vcov") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());

  REQUIRE(r->psi.rows() == x.rows());
  REQUIRE(r->psi.cols() == 3);
  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - r->vcov).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("estimate_fiml_gwise: complete data matches complete g-wise estimator") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance.has_value());

  auto complete = ms::estimate_gwise(x, *distance);
  auto fiml = ms::estimate_fiml_gwise(x, *distance, EmOptions{});
  REQUIRE(complete.has_value());
  REQUIRE(fiml.has_value());

  CHECK((complete->estimates - fiml->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((complete->vcov - fiml->vcov).cwiseAbs().maxCoeff() < 1e-8);
  REQUIRE(fiml->psi.rows() == x.rows());
  REQUIRE(fiml->psi.cols() == 2);
  const RealMat psi_vcov =
      (fiml->psi.transpose() * fiml->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - fiml->vcov).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("estimate_fiml_gwise: g=2 nominal matches categorical FIML") {
  IntMat x = twelve_subject_3rater_3cat();
  auto distance = ms::loss::frechet_nominal_distance(3);
  auto weights = ms::loss::identity_weights(3);
  REQUIRE(distance.has_value());
  REQUIRE(weights.has_value());

  auto gwise = ms::estimate_fiml_gwise(x, *distance, EmOptions{}, ms::GwiseOptions{2});
  auto pairwise = ms::estimate_fiml(x, *weights, EmOptions{});
  REQUIRE(gwise.has_value());
  REQUIRE(pairwise.has_value());

  CHECK(std::abs(gwise->estimates(0) - pairwise->estimates(0)) < 1e-8);
  CHECK(std::abs(gwise->estimates(1) - pairwise->estimates(1)) < 1e-8);
  CHECK((gwise->vcov - pairwise->vcov.block(0, 0, 2, 2)).cwiseAbs().maxCoeff() < 1e-6);
}

TEST_CASE("estimate_fiml_gwise: variance is symmetric PSD and psi reconstructs it") {
  IntMat x = twelve_subject_3rater_3cat();
  auto distance = ms::loss::hubert_categorical_distance(3);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_fiml_gwise(x, *distance, EmOptions{});
  REQUIRE(r.has_value());
  const RealMat& V = r->vcov;
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-10);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-8);

  REQUIRE(r->psi.rows() == x.rows());
  REQUIRE(r->psi.cols() == 2);
  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - r->vcov).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("estimate_fiml: info_rcond affects Louis variance, not estimates") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);

  EmOptions keep{};
  keep.info_rcond = 0.0;
  auto r_keep = ms::estimate_fiml(x, *W, keep);
  REQUIRE(r_keep.has_value());

  EmOptions drop{};
  drop.info_rcond = 1.0;
  auto r_drop = ms::estimate_fiml(x, *W, drop);
  REQUIRE(r_drop.has_value());

  CHECK((r_keep->estimates - r_drop->estimates).cwiseAbs().maxCoeff() < 1e-12);
  CHECK(r_drop->vcov.norm() < r_keep->vcov.norm());
}

TEST_CASE("estimate_fiml: perfect agreement converges to kappa = 1") {
  IntMat x(5, 2);
  x << 0, 0,
       1, 1,
       0, 0,
       1, 1,
       0, 0;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < 1e-6);
  CHECK(std::abs(r->estimates(1) - 1.0) < 1e-6);
  CHECK(std::abs(r->estimates(2) - 1.0) < 1e-6);
}

TEST_CASE("estimate_fiml_gwise: perfect agreement converges to kappa = 1") {
  IntMat x(5, 3);
  x << 0, 0, 0,
       1, 1, 1,
       0, 0, 0,
       1, 1, 1,
       0, 0, 0;
  auto distance = ms::loss::frechet_nominal_distance(2);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_fiml_gwise(x, *distance, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < 1e-6);
  CHECK(std::abs(r->estimates(1) - 1.0) < 1e-6);
}

TEST_CASE("estimate_fiml: too few raters -> invalid_argument") {
  IntMat x(3, 1);
  x << 0, 1, 0;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}

TEST_CASE("estimate_fiml_gwise: rejects invalid arguments") {
  auto distance = ms::loss::frechet_nominal_distance(3);
  REQUIRE(distance.has_value());

  IntMat one_rater(3, 1);
  one_rater << 0, 1, 0;
  auto too_few = ms::estimate_fiml_gwise(one_rater, *distance, EmOptions{});
  REQUIRE(!too_few.has_value());
  CHECK(too_few.error() == ms::Error::invalid_argument);

  IntMat out_of_range = twelve_subject_3rater_3cat();
  out_of_range(0, 0) = 3;
  auto bad_category = ms::estimate_fiml_gwise(out_of_range, *distance, EmOptions{});
  REQUIRE(!bad_category.has_value());
  CHECK(bad_category.error() == ms::Error::invalid_argument);

  auto distance5 = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance5.has_value());
  auto oversized =
      ms::estimate_fiml_gwise(
          frechet_fixture_int(), *distance5, EmOptions{}, ms::GwiseOptions{5, 100});
  REQUIRE(!oversized.has_value());
  CHECK(oversized.error() == ms::Error::not_supported);
}

TEST_CASE("estimate_fiml: out-of-range category -> invalid_argument") {
  IntMat x(3, 2);
  x << 0, 1,
       5, 0,
       0, 1;
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::invalid_argument);
}
