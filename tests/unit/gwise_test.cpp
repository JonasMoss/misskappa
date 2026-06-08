#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

using misskappa::IntMat;
using misskappa::RealMat;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

IntMat frechet_fixture_int() {
  IntMat x(4, 5);
  x << 0, 0, 1, 0, 0,
       0, 1, 2, 1, 1,
       1, 0, 0, 0, 0,
       1, 2, 3, 3, 4;
  return x;
}

RealMat frechet_fixture_real() {
  return frechet_fixture_int().cast<double>() + RealMat::Constant(4, 5, 1.0);
}

RealMat missing_3rater_real_fixture() {
  const double na = std::numeric_limits<double>::quiet_NaN();
  RealMat x(12, 3);
  x << 0.0, 0.0, 0.0,
       1.0, 1.0, na,
       2.0, 2.0, 2.0,
       0.0, 0.0, 1.0,
       na, 2.0, 1.0,
       2.0, 1.0, 2.0,
       0.0, 1.0, 0.0,
       1.0, 1.0, 2.0,
       2.0, 2.0, 2.0,
       0.0, 0.0, 0.0,
       1.0, 0.0, 1.0,
       2.0, 2.0, 1.0;
  return x;
}

IntMat missing_3rater_fixture() {
  IntMat x(12, 3);
  x <<
    0, 0, 0,
    1, 1, ms::na_code,
    2, 2, 2,
    0, 0, 1,
    ms::na_code, 2, 1,
    2, 1, 2,
    0, 1, 0,
    1, 1, 2,
    2, 2, 2,
    0, 0, 0,
    1, 0, 1,
    2, 2, 1;
  return x;
}

void check_symmetric_psd(const RealMat& V) {
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-12);
  Eigen::SelfAdjointEigenSolver<RealMat> es(V);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-10);
}

}  // namespace

TEST_CASE("G-wise nominal Frechet: g=2 matches pairwise identity kappa") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  auto weights = ms::loss::identity_weights(5);
  REQUIRE(distance.has_value());
  REQUIRE(weights.has_value());

  auto gwise = ms::estimate_gwise(x, *distance, ms::GwiseOptions{2});
  auto pairwise = ms::estimate_available(x, *weights);
  REQUIRE(gwise.has_value());
  REQUIRE(pairwise.has_value());

  CHECK(std::abs(gwise->estimates(0) - pairwise->estimates(0)) < tol);
  CHECK(std::abs(gwise->estimates(1) - pairwise->estimates(1)) < tol);
  CHECK((gwise->vcov - pairwise->vcov.block(0, 0, 2, 2)).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("G-wise nominal Frechet IPW: complete data matches complete estimator") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance.has_value());

  auto complete = ms::estimate_gwise(x, *distance);
  auto ipw = ms::estimate_ipw_gwise(x, *distance);
  REQUIRE(complete.has_value());
  REQUIRE(ipw.has_value());

  CHECK((complete->estimates - ipw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((complete->vcov - ipw->vcov).cwiseAbs().maxCoeff() < 1e-12);
  CHECK((complete->psi - ipw->psi).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("G-wise nominal Frechet IPW: g=2 matches categorical IPW") {
  IntMat x = missing_3rater_fixture();
  auto distance = ms::loss::frechet_nominal_distance(3);
  auto weights = ms::loss::identity_weights(3);
  REQUIRE(distance.has_value());
  REQUIRE(weights.has_value());

  auto gwise = ms::estimate_ipw_gwise(x, *distance, ms::GwiseOptions{2});
  auto pairwise = ms::estimate_ipw(x, *weights);
  REQUIRE(gwise.has_value());
  REQUIRE(pairwise.has_value());

  CHECK(std::abs(gwise->estimates(0) - pairwise->estimates(0)) < tol);
  CHECK(std::abs(gwise->estimates(1) - pairwise->estimates(1)) < tol);
  CHECK((gwise->vcov - pairwise->vcov.block(0, 0, 2, 2)).cwiseAbs().maxCoeff() < 1e-6);

  const RealMat psi_vcov =
      (gwise->psi.transpose() * gwise->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - gwise->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("G-wise nominal Frechet: all-rater mode disagreement fixture") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise(x, *distance);
  REQUIRE(r.has_value());
  REQUIRE(r->estimates.size() == 2);
  CHECK(std::abs(r->estimates(0) - 0.23353293413173815) < tol);
  CHECK(std::abs(r->estimates(1) - 0.21735788407113500) < tol);
  CHECK(r->psi.rows() == x.rows());
  CHECK(r->psi.cols() == 2);
}

TEST_CASE("G-wise absolute Frechet: all-rater median disagreement fixture") {
  RealMat x = frechet_fixture_real();
  auto distance = ms::loss::frechet_absolute_distance();
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise_continuous(x, *distance);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 0.45877378435517835) < tol);
  CHECK(std::abs(r->estimates(1) - 0.44821878125323300) < tol);

  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - r->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("G-wise absolute Frechet IPW: complete data matches complete estimator") {
  RealMat x = frechet_fixture_real();
  auto distance = ms::loss::frechet_absolute_distance();
  REQUIRE(distance.has_value());

  auto complete = ms::estimate_gwise_continuous(x, *distance);
  auto ipw = ms::estimate_ipw_gwise_continuous(x, *distance);
  REQUIRE(complete.has_value());
  REQUIRE(ipw.has_value());

  CHECK((complete->estimates - ipw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((complete->vcov - ipw->vcov).cwiseAbs().maxCoeff() < 1e-12);
  CHECK((complete->psi - ipw->psi).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("G-wise continuous IPW: g=2 matches pairwise continuous IPW") {
  RealMat x = missing_3rater_real_fixture();

  auto absolute_distance = ms::loss::frechet_absolute_distance();
  auto absolute_loss = ms::loss::linear_loss(0.0, 2.0);
  REQUIRE(absolute_distance.has_value());
  REQUIRE(absolute_loss.has_value());
  auto gwise_absolute =
      ms::estimate_ipw_gwise_continuous(x, *absolute_distance, ms::GwiseOptions{2});
  auto pairwise_absolute = ms::estimate_ipw_continuous(x, *absolute_loss);
  REQUIRE(gwise_absolute.has_value());
  REQUIRE(pairwise_absolute.has_value());
  CHECK((gwise_absolute->estimates - pairwise_absolute->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((gwise_absolute->vcov - pairwise_absolute->vcov).cwiseAbs().maxCoeff() < 1e-10);

  auto quadratic_distance = ms::loss::frechet_quadratic_distance();
  auto quadratic_loss = ms::loss::quadratic_loss(0.0, 2.0);
  REQUIRE(quadratic_distance.has_value());
  REQUIRE(quadratic_loss.has_value());
  auto gwise_quadratic =
      ms::estimate_ipw_gwise_continuous(x, *quadratic_distance, ms::GwiseOptions{2});
  auto pairwise_quadratic = ms::estimate_ipw_continuous(x, *quadratic_loss);
  REQUIRE(gwise_quadratic.has_value());
  REQUIRE(pairwise_quadratic.has_value());
  CHECK((gwise_quadratic->estimates - pairwise_quadratic->estimates).cwiseAbs().maxCoeff() < tol);

  auto hubert_distance = ms::loss::hubert_continuous_distance();
  auto identity_loss = ms::loss::identity_loss();
  REQUIRE(hubert_distance.has_value());
  REQUIRE(identity_loss.has_value());
  auto gwise_hubert =
      ms::estimate_ipw_gwise_continuous(x, *hubert_distance, ms::GwiseOptions{2});
  auto pairwise_identity = ms::estimate_ipw_continuous(x, *identity_loss);
  REQUIRE(gwise_hubert.has_value());
  REQUIRE(pairwise_identity.has_value());
  CHECK((gwise_hubert->estimates - pairwise_identity->estimates).cwiseAbs().maxCoeff() < tol);

  const RealMat psi_vcov = (gwise_absolute->psi.transpose() * gwise_absolute->psi)
                           / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - gwise_absolute->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("G-wise Hubert disagreement: perfect agreement returns one") {
  IntMat x(3, 3);
  x << 0, 0, 0,
       1, 1, 1,
       2, 2, 2;
  auto distance = ms::loss::hubert_categorical_distance(3);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise(x, *distance);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
}

TEST_CASE("G-wise continuous IPW: perfect observed agreement returns one") {
  const double na = std::numeric_limits<double>::quiet_NaN();
  RealMat x(4, 3);
  x << 0.0, 0.0, 0.0,
       1.0, 1.0, na,
       2.0, 2.0, 2.0,
       na, 1.0, 1.0;
  auto distance = ms::loss::frechet_absolute_distance();
  REQUIRE(distance.has_value());

  auto r = ms::estimate_ipw_gwise_continuous(x, *distance);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
}

TEST_CASE("G-wise nominal Frechet IPW: perfect observed agreement returns one") {
  IntMat x(4, 3);
  x << 0, 0, 0,
       1, 1, ms::na_code,
       2, 2, 2,
       ms::na_code, 1, 1;
  auto distance = ms::loss::frechet_nominal_distance(3);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_ipw_gwise(x, *distance);
  REQUIRE(r.has_value());
  CHECK(std::abs(r->estimates(0) - 1.0) < tol);
  CHECK(std::abs(r->estimates(1) - 1.0) < tol);
}

TEST_CASE("G-wise estimator rejects missing and oversized direct jobs") {
  IntMat x = frechet_fixture_int();
  auto distance = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance.has_value());

  x(0, 0) = ms::na_code;
  auto missing = ms::estimate_gwise(x, *distance);
  REQUIRE(!missing.has_value());
  CHECK(missing.error() == ms::Error::invalid_argument);

  x = frechet_fixture_int();
  auto oversized = ms::estimate_gwise(x, *distance, ms::GwiseOptions{5, 100});
  REQUIRE(!oversized.has_value());
  CHECK(oversized.error() == ms::Error::not_supported);
}

TEST_CASE("G-wise nominal Frechet IPW: rejects invalid or unidentified missingness patterns") {
  auto distance = ms::loss::frechet_nominal_distance(3);
  REQUIRE(distance.has_value());

  IntMat out_of_range = missing_3rater_fixture();
  out_of_range(0, 0) = 3;
  auto bad_category = ms::estimate_ipw_gwise(out_of_range, *distance);
  REQUIRE(!bad_category.has_value());
  CHECK(bad_category.error() == ms::Error::invalid_argument);

  IntMat zero_rater(3, 3);
  zero_rater << 0, 0, ms::na_code,
                1, 1, ms::na_code,
                2, 2, ms::na_code;
  auto singular_rater = ms::estimate_ipw_gwise(zero_rater, *distance, ms::GwiseOptions{2});
  REQUIRE(!singular_rater.has_value());
  CHECK(singular_rater.error() == ms::Error::singular_weight);

  IntMat no_joint(3, 3);
  no_joint << 0, ms::na_code, ms::na_code,
              ms::na_code, 1, ms::na_code,
              ms::na_code, ms::na_code, 2;
  auto no_tuple = ms::estimate_ipw_gwise(no_joint, *distance);
  REQUIRE(!no_tuple.has_value());
  CHECK(no_tuple.error() == ms::Error::singular_weight);

  auto distance5 = ms::loss::frechet_nominal_distance(5);
  REQUIRE(distance5.has_value());
  auto oversized =
      ms::estimate_ipw_gwise(frechet_fixture_int(), *distance5, ms::GwiseOptions{5, 100});
  REQUIRE(!oversized.has_value());
  CHECK(oversized.error() == ms::Error::not_supported);
}

TEST_CASE("G-wise continuous IPW: rejects unidentified missingness patterns and oversized jobs") {
  const double na = std::numeric_limits<double>::quiet_NaN();
  auto distance = ms::loss::frechet_absolute_distance();
  REQUIRE(distance.has_value());

  RealMat zero_rater(3, 3);
  zero_rater << 0.0, 0.0, na,
                1.0, 1.0, na,
                2.0, 2.0, na;
  auto singular_rater =
      ms::estimate_ipw_gwise_continuous(zero_rater, *distance, ms::GwiseOptions{2});
  REQUIRE(!singular_rater.has_value());
  CHECK(singular_rater.error() == ms::Error::singular_weight);

  RealMat no_joint(3, 3);
  no_joint << 0.0, na, na,
              na, 1.0, na,
              na, na, 2.0;
  auto no_tuple = ms::estimate_ipw_gwise_continuous(no_joint, *distance);
  REQUIRE(!no_tuple.has_value());
  CHECK(no_tuple.error() == ms::Error::singular_weight);

  auto oversized =
      ms::estimate_ipw_gwise_continuous(
          frechet_fixture_real(), *distance, ms::GwiseOptions{5, 100});
  REQUIRE(!oversized.has_value());
  CHECK(oversized.error() == ms::Error::not_supported);
}

TEST_CASE("G-wise categorical chance uses finite support instead of n^g item tuples") {
  IntMat x(200, 3);
  for (int i = 0; i < static_cast<int>(x.rows()); ++i) {
    x(i, 0) = i % 2;
    x(i, 1) = (i / 2) % 2;
    x(i, 2) = (i + 1) % 2;
  }
  auto distance = ms::loss::frechet_nominal_distance(2);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise(x, *distance, ms::GwiseOptions{3, 100});
  REQUIRE(r.has_value());
  CHECK(std::isfinite(r->estimates(0)));
  CHECK(std::isfinite(r->estimates(1)));
}

TEST_CASE("G-wise covariance is symmetric and PSD") {
  RealMat x = frechet_fixture_real();
  auto distance = ms::loss::frechet_quadratic_distance();
  REQUIRE(distance.has_value());

  auto r = ms::estimate_gwise_continuous(x, *distance);
  REQUIRE(r.has_value());
  CHECK((r->vcov - r->vcov.transpose()).cwiseAbs().maxCoeff() < 1e-12);

  Eigen::SelfAdjointEigenSolver<RealMat> es(r->vcov);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-10);
}

TEST_CASE("G-wise nominal Frechet IPW covariance is symmetric and PSD") {
  IntMat x = missing_3rater_fixture();
  auto distance = ms::loss::frechet_nominal_distance(3);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_ipw_gwise(x, *distance);
  REQUIRE(r.has_value());
  check_symmetric_psd(r->vcov);
}

TEST_CASE("G-wise continuous IPW covariance is symmetric and PSD") {
  RealMat x = missing_3rater_real_fixture();
  auto distance = ms::loss::frechet_absolute_distance();
  REQUIRE(distance.has_value());

  auto r = ms::estimate_ipw_gwise_continuous(x, *distance);
  REQUIRE(r.has_value());
  check_symmetric_psd(r->vcov);
}
