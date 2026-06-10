#include "doctest.h"

#include "misskappa/estimate.hpp"

#include <cmath>
#include <limits>

using misskappa::EmOptions;
using misskappa::IntMat;
using misskappa::na_code;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

IntMat complete_items() {
  IntMat x(6, 3);
  x << 0, 0, 1,
       1, 1, 1,
       2, 1, 2,
       1, 2, 2,
       0, 1, 0,
       2, 2, 1;
  return x;
}

IntMat missing_items() {
  IntMat x(7, 3);
  x << 0, 0, 1,
       1, na_code, 1,
       2, 1, 2,
       na_code, 2, 2,
       0, 1, 0,
       2, 2, na_code,
       1, 1, 1;
  return x;
}

IntMat identified_missing_items() {
  IntMat x(36, 3);
  Eigen::Index row = 0;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      for (int c = 0; c < 3; ++c) {
        x(row, 0) = a;
        x(row, 1) = b;
        x(row, 2) = c;
        ++row;
      }
    }
  }
  for (int a = 0; a < 3; ++a) {
    x(row, 0) = a;
    x(row, 1) = na_code;
    x(row, 2) = (a + 1) % 3;
    ++row;
    x(row, 0) = na_code;
    x(row, 1) = a;
    x(row, 2) = (a + 2) % 3;
    ++row;
    x(row, 0) = a;
    x(row, 1) = (a + 1) % 3;
    x(row, 2) = na_code;
    ++row;
  }
  return x;
}

RealVec three_scores() {
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  return v;
}

RealMat continuous_items() {
  RealMat x(7, 3);
  x << 1.0, 1.2, 1.1,
       2.0, 2.1, 2.0,
       3.0, 2.8, 3.2,
       2.0, 2.2, 2.4,
       1.0, 1.5, 1.2,
       3.0, 2.9, 2.7,
       2.0, 2.0, 2.1;
  return x;
}

RealMat missing_continuous_items() {
  RealMat x = continuous_items();
  x(1, 1) = std::numeric_limits<double>::quiet_NaN();
  x(3, 0) = std::numeric_limits<double>::quiet_NaN();
  x(5, 2) = std::numeric_limits<double>::quiet_NaN();
  return x;
}

double complete_alpha_by_covariance(const IntMat& x, const RealVec& values) {
  const int n = static_cast<int>(x.rows());
  const int R = static_cast<int>(x.cols());
  RealMat y(n, R);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < R; ++j) y(i, j) = values(x(i, j));
  }
  RealVec mu = y.colwise().mean().transpose();
  RealMat sigma = RealMat::Zero(R, R);
  for (int i = 0; i < n; ++i) {
    const RealVec centered = y.row(i).transpose() - mu;
    sigma.noalias() += centered * centered.transpose();
  }
  sigma /= static_cast<double>(n);
  const double t1 = sigma.sum();
  const double t2 = sigma.diagonal().sum();
  return (static_cast<double>(R) / static_cast<double>(R - 1)) * (1.0 - t2 / t1);
}

double complete_alpha_by_covariance(const RealMat& x) {
  const int n = static_cast<int>(x.rows());
  const int R = static_cast<int>(x.cols());
  RealVec mu = x.colwise().mean().transpose();
  RealMat sigma = RealMat::Zero(R, R);
  for (int i = 0; i < n; ++i) {
    const RealVec centered = x.row(i).transpose() - mu;
    sigma.noalias() += centered * centered.transpose();
  }
  sigma /= static_cast<double>(n);
  const double t1 = sigma.sum();
  const double t2 = sigma.diagonal().sum();
  return (static_cast<double>(R) / static_cast<double>(R - 1)) * (1.0 - t2 / t1);
}

void check_influence_reconstructs_vcov(const ms::Estimation& e, int n) {
  REQUIRE(e.psi.rows() == n);
  REQUIRE(e.psi.cols() == e.estimates.size());
  const RealMat psi_vcov =
      (e.psi.transpose() * e.psi) / std::pow(static_cast<double>(n), 2);
  CHECK((psi_vcov - e.vcov).cwiseAbs().maxCoeff() < 1e-10);
}

}  // namespace

TEST_CASE("estimate_alpha_available: complete data matches covariance alpha") {
  const IntMat x = complete_items();
  const RealVec values = three_scores();
  auto r = ms::estimate_alpha_available(x, values);
  REQUIRE(r.has_value());

  CHECK(std::abs(r->estimates(0) - complete_alpha_by_covariance(x, values)) < 1e-12);
  CHECK(r->vcov.rows() == 1);
  CHECK(r->vcov.cols() == 1);
  CHECK(r->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*r, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_available_continuous: complete data matches covariance alpha") {
  const RealMat x = continuous_items();
  auto r = ms::estimate_alpha_available_continuous(x);
  REQUIRE(r.has_value());

  CHECK(std::abs(r->estimates(0) - complete_alpha_by_covariance(x)) < 1e-12);
  CHECK(r->vcov.rows() == 1);
  CHECK(r->vcov.cols() == 1);
  CHECK(r->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*r, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_available_continuous: missing fixture is finite") {
  const RealMat x = missing_continuous_items();
  auto r = ms::estimate_alpha_available_continuous(x);
  REQUIRE(r.has_value());

  CHECK(std::isfinite(r->estimates(0)));
  CHECK(r->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*r, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_available: scored categorical path matches continuous path") {
  const IntMat x = complete_items();
  const RealVec values = three_scores();
  RealMat y(x.rows(), x.cols());
  for (int i = 0; i < x.rows(); ++i) {
    for (int j = 0; j < x.cols(); ++j) y(i, j) = values(x(i, j));
  }

  auto cat = ms::estimate_alpha_available(x, values);
  auto con = ms::estimate_alpha_available_continuous(y);
  REQUIRE(cat.has_value());
  REQUIRE(con.has_value());

  CHECK(std::abs(cat->estimates(0) - con->estimates(0)) < 1e-12);
  CHECK((cat->vcov - con->vcov).cwiseAbs().maxCoeff() < 1e-12);
  CHECK((cat->psi - con->psi).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("estimate_alpha_fiml: complete data matches available alpha") {
  const IntMat x = complete_items();
  const RealVec values = three_scores();
  auto av = ms::estimate_alpha_available(x, values);
  REQUIRE(av.has_value());
  auto ml = ms::estimate_alpha_fiml(x, values, EmOptions{});
  REQUIRE(ml.has_value());

  CHECK(std::abs(ml->estimates(0) - av->estimates(0)) < 1e-10);
  CHECK(ml->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*ml, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_fiml: sparse missing categorical fixture reports non-identification") {
  const IntMat x = missing_items();
  const RealVec values = three_scores();
  auto av = ms::estimate_alpha_available(x, values);
  REQUIRE(av.has_value());
  auto ml = ms::estimate_alpha_fiml(x, values, EmOptions{});
  REQUIRE(!ml.has_value());

  CHECK(std::isfinite(av->estimates(0)));
  CHECK(av->vcov(0, 0) >= -1e-12);
  CHECK(ml.error() == ms::Error::not_identified);
  check_influence_reconstructs_vcov(*av, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha_fiml: identified missing categorical fixture is finite") {
  const IntMat x = identified_missing_items();
  const RealVec values = three_scores();
  auto ml = ms::estimate_alpha_fiml(x, values, EmOptions{});
  REQUIRE(ml.has_value());

  CHECK(std::isfinite(ml->estimates(0)));
  CHECK(ml->vcov(0, 0) >= -1e-12);
  check_influence_reconstructs_vcov(*ml, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_alpha: invalid inputs are rejected") {
  const RealVec values = three_scores();

  IntMat too_few_items(3, 1);
  too_few_items << 0, 1, 2;
  auto few = ms::estimate_alpha_available(too_few_items, values);
  REQUIRE(!few.has_value());
  CHECK(few.error() == ms::Error::invalid_argument);

  IntMat bad = complete_items();
  bad(0, 0) = 5;
  auto out_of_range = ms::estimate_alpha_fiml(bad, values, EmOptions{});
  REQUIRE(!out_of_range.has_value());
  CHECK(out_of_range.error() == ms::Error::invalid_argument);
}

TEST_CASE("estimate_alpha: incomplete item co-observation graph is not identified") {
  const RealVec values = three_scores();
  IntMat categorical(6, 3);
  categorical << 0, 0, na_code,
                 0, 1, na_code,
                 1, 0, na_code,
                 na_code, 0, 0,
                 na_code, 0, 1,
                 na_code, 1, 1;
  RealMat continuous(6, 3);
  continuous << 1.0, 1.0, std::numeric_limits<double>::quiet_NaN(),
                1.0, 2.0, std::numeric_limits<double>::quiet_NaN(),
                2.0, 1.0, std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0,
                std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0,
                std::numeric_limits<double>::quiet_NaN(), 2.0, 2.0;

  auto cat_av = ms::estimate_alpha_available(categorical, values);
  auto cat_ml = ms::estimate_alpha_fiml(categorical, values, EmOptions{});
  auto con_av = ms::estimate_alpha_available_continuous(continuous);
  REQUIRE(!cat_av.has_value());
  REQUIRE(!cat_ml.has_value());
  REQUIRE(!con_av.has_value());
  CHECK(cat_av.error() == ms::Error::not_identified);
  CHECK(cat_ml.error() == ms::Error::not_identified);
  CHECK(con_av.error() == ms::Error::not_identified);
}
