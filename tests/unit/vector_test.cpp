#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <cmath>
#include <limits>

using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;
const double na_d = std::numeric_limits<double>::quiet_NaN();

RealMat vector_fixture_complete() {
  // n = 5, R = 3, p = 2, rater-major columns.
  RealMat x(5, 6);
  x <<
    0, 1,  0, 1,  1, 1,
    1, 0,  1, 0,  1, 1,
    0, 0,  1, 0,  0, 0,
    1, 1,  1, 0,  0, 1,
    0, 1,  0, 0,  0, 1;
  return x;
}

RealVec two_feature_weights() {
  RealVec w(2);
  w << 1.0, 2.0;
  return w;
}

RealVec squared_moment_kappa(const RealMat& x, int R, int features, const RealVec& w) {
  const int n = static_cast<int>(x.rows());
  RealVec mu = x.colwise().mean().transpose();
  RealMat centered = x.rowwise() - mu.transpose();
  RealMat Sigma = (centered.transpose() * centered) / static_cast<double>(n);

  double T = 0.0;
  double B = 0.0;
  double G = 0.0;
  for (int l = 0; l < features; ++l) {
    double mu_bar = 0.0;
    for (int r = 0; r < R; ++r) mu_bar += mu(r * features + l);
    mu_bar /= static_cast<double>(R);

    for (int r = 0; r < R; ++r) {
      const int cr = r * features + l;
      T += w(l) * Sigma(cr, cr);
      const double delta = mu(cr) - mu_bar;
      G += w(l) * delta * delta;
      for (int s = 0; s < R; ++s) {
        const int cs = s * features + l;
        B += w(l) * Sigma(cr, cs);
      }
    }
  }

  RealVec out(2);
  out(0) = (B - T) / ((R - 1.0) * T + R * G);
  out(1) = (B - T - G) / ((R - 1.0) * (T + G));
  return out;
}

}  // namespace

TEST_CASE("Vector hamming with one feature matches scalar identity loss") {
  RealMat x(6, 3);
  x <<
    0, 0, 1,
    1, 1, 1,
    0, 1, 0,
    1, na_d, 0,
    0, 0, na_d,
    1, 0, 1;

  RealVec w(1);
  w << 1.0;
  auto vloss = ms::loss::hamming_vector_loss(w);
  auto closs = ms::loss::identity_loss();
  REQUIRE(vloss.has_value());
  REQUIRE(closs.has_value());

  auto vec = ms::estimate_pairwise_vector(x, 1, *vloss);
  auto con = ms::estimate_available_continuous(x, *closs);
  REQUIRE(vec.has_value());
  REQUIRE(con.has_value());

  CHECK((vec->estimates - con->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((vec->vcov - con->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("Vector complete data: pairwise and IPW collapse") {
  RealMat x = vector_fixture_complete();
  auto loss = ms::loss::absolute_vector_loss(two_feature_weights());
  REQUIRE(loss.has_value());

  auto pairwise = ms::estimate_pairwise_vector(x, 2, *loss);
  auto ipw = ms::estimate_ipw_vector(x, 2, *loss);
  REQUIRE(pairwise.has_value());
  REQUIRE(ipw.has_value());

  CHECK((pairwise->estimates - ipw->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK((pairwise->vcov - ipw->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("Vector squared loss matches diagonal mean-covariance contraction") {
  RealMat x = vector_fixture_complete();
  RealVec w = two_feature_weights();
  auto loss = ms::loss::squared_vector_loss(w);
  REQUIRE(loss.has_value());

  auto fit = ms::estimate_pairwise_vector(x, 2, *loss);
  REQUIRE(fit.has_value());
  RealVec expected = squared_moment_kappa(x, 3, 2, w);

  CHECK(std::abs(fit->estimates(0) - expected(0)) < tol);
  CHECK(std::abs(fit->estimates(1) - expected(1)) < tol);
}

TEST_CASE("Vector component-wise missing fixture has IF covariance") {
  RealMat x = vector_fixture_complete();
  x(0, 1) = na_d;
  x(1, 4) = na_d;
  x(3, 2) = na_d;
  auto loss = ms::loss::rms_vector_loss(two_feature_weights());
  REQUIRE(loss.has_value());

  auto fit = ms::estimate_ipw_vector(x, 2, *loss);
  REQUIRE(fit.has_value());
  REQUIRE(fit->estimates.size() == 2);
  CHECK(std::isfinite(fit->estimates(0)));
  CHECK(std::isfinite(fit->estimates(1)));
  REQUIRE(fit->psi.rows() == x.rows());
  REQUIRE(fit->psi.cols() == 2);

  const RealMat psi_vcov =
      (fit->psi.transpose() * fit->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - fit->vcov).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("Vector estimator rejects invalid inputs") {
  RealMat x = vector_fixture_complete();
  auto good_loss = ms::loss::hamming_vector_loss(two_feature_weights());
  REQUIRE(good_loss.has_value());

  auto bad_features = ms::estimate_pairwise_vector(x, 0, *good_loss);
  REQUIRE(!bad_features.has_value());
  CHECK(bad_features.error() == ms::Error::invalid_argument);

  auto bad_cols = ms::estimate_pairwise_vector(x, 4, *good_loss);
  REQUIRE(!bad_cols.has_value());
  CHECK(bad_cols.error() == ms::Error::dimension_mismatch);

  RealVec zero_w(2);
  zero_w << 0.0, 0.0;
  auto bad_loss = ms::loss::squared_vector_loss(zero_w);
  REQUIRE(!bad_loss.has_value());
  CHECK(bad_loss.error() == ms::Error::invalid_argument);

  x.col(1).array() = na_d;
  auto singular = ms::estimate_ipw_vector(x, 2, *good_loss);
  REQUIRE(!singular.has_value());
  CHECK(singular.error() == ms::Error::singular_weight);
}
