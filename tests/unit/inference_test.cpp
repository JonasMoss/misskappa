#include "doctest.h"

#include "misskappa/estimate.hpp"
#include "misskappa/inference.hpp"
#include "misskappa/loss.hpp"

#include <cmath>

namespace ms = misskappa;

namespace {

constexpr double tol = 1e-11;

ms::IntMat delete_int_row(const ms::IntMat& x, int omitted) {
  ms::IntMat out(x.rows() - 1, x.cols());
  int target = 0;
  for (int i = 0; i < x.rows(); ++i) {
    if (i == omitted) continue;
    out.row(target++) = x.row(i);
  }
  return out;
}

ms::RealMat delete_real_row(const ms::RealMat& x, int omitted) {
  ms::RealMat out(x.rows() - 1, x.cols());
  int target = 0;
  for (int i = 0; i < x.rows(); ++i) {
    if (i == omitted) continue;
    out.row(target++) = x.row(i);
  }
  return out;
}

}  // namespace

TEST_CASE("categorical complete-data delete-one matches exact refits") {
  ms::IntMat x(7, 3);
  x << 0, 0, 1,
       1, 1, 1,
       2, 1, 2,
       1, 2, 2,
       0, 1, 0,
       2, 2, 1,
       0, 2, 1;
  ms::RealVec values(3);
  values << 1.0, 2.0, 4.0;
  auto weights = ms::loss::linear_weights(3, values);
  REQUIRE(weights.has_value());

  auto summary = ms::delete_one_kappa_complete(x, *weights);
  REQUIRE(summary.has_value());
  auto full = ms::estimate_available(x, *weights);
  REQUIRE(full.has_value());
  CHECK((summary->coefficient.estimate - full->estimates).cwiseAbs().maxCoeff() < tol);
  REQUIRE(summary->item_disagreement.size() == x.rows());

  for (int i = 0; i < x.rows(); ++i) {
    const ms::IntMat reduced = delete_int_row(x, i);
    auto fit = ms::estimate_available(reduced, *weights);
    REQUIRE(fit.has_value());
    CHECK((summary->coefficient.leave_one_out.row(i).transpose() - fit->estimates)
              .cwiseAbs().maxCoeff() < tol);
  }
}

TEST_CASE("continuous complete-data delete-one matches exact refits") {
  ms::RealMat x(6, 3);
  x << 1.0, 1.2, 1.1,
       2.0, 2.1, 2.0,
       3.0, 2.8, 3.2,
       2.0, 2.2, 2.4,
       1.0, 1.5, 1.2,
       3.0, 2.9, 2.7;
  auto loss = ms::loss::quadratic_loss(1.0, 3.2);
  REQUIRE(loss.has_value());

  auto summary = ms::delete_one_kappa_complete(x, *loss);
  REQUIRE(summary.has_value());
  auto full = ms::estimate_available_continuous(x, *loss);
  REQUIRE(full.has_value());
  CHECK((summary->coefficient.estimate - full->estimates).cwiseAbs().maxCoeff() < tol);

  for (int i = 0; i < x.rows(); ++i) {
    const ms::RealMat reduced = delete_real_row(x, i);
    auto fit = ms::estimate_available_continuous(reduced, *loss);
    REQUIRE(fit.has_value());
    CHECK((summary->coefficient.leave_one_out.row(i).transpose() - fit->estimates)
              .cwiseAbs().maxCoeff() < tol);
  }
}

TEST_CASE("count-data delete-one matches exact Fleiss-Cuzick refits") {
  ms::IntMat x(7, 3);
  x << 4, 0, 0,
       0, 4, 0,
       0, 0, 4,
       3, 1, 0,
       1, 3, 0,
       0, 3, 1,
       2, 1, 1;
  ms::RealVec values(3);
  values << 1.0, 2.0, 3.0;
  auto weights = ms::loss::quadratic_weights(3, values);
  REQUIRE(weights.has_value());

  auto summary = ms::delete_one_kappa_counts_complete(x, *weights);
  REQUIRE(summary.has_value());
  auto full = ms::estimate_available_counts(x, *weights);
  REQUIRE(full.has_value());
  CHECK((summary->coefficient.estimate - full->estimates).cwiseAbs().maxCoeff() < tol);
  REQUIRE(summary->item_disagreement.size() == x.rows());

  for (int i = 0; i < x.rows(); ++i) {
    const ms::IntMat reduced = delete_int_row(x, i);
    auto fit = ms::estimate_available_counts(reduced, *weights);
    REQUIRE(fit.has_value());
    CHECK((summary->coefficient.leave_one_out.row(i).transpose() - fit->estimates)
              .cwiseAbs().maxCoeff() < tol);
  }
}

TEST_CASE("unit-weighted count delete-one matches exact variable-row refits") {
  ms::IntMat x(5, 3);
  x << 5, 0, 0,
       1, 1, 0,
       0, 2, 1,
       1, 1, 2,
       0, 0, 3;
  auto weights = ms::loss::identity_weights(3);
  REQUIRE(weights.has_value());

  auto summary = ms::delete_one_kappa_counts_complete(
      x, *weights, ms::CountWeighting::unit_weighted);
  REQUIRE(summary.has_value());
  auto full = ms::estimate_counts(
      x, *weights, ms::CountWeighting::unit_weighted);
  REQUIRE(full.has_value());
  CHECK((summary->coefficient.estimate - full->estimates).cwiseAbs().maxCoeff() < tol);
  CHECK(summary->item_disagreement.size() == 0);

  for (int i = 0; i < x.rows(); ++i) {
    const ms::IntMat reduced = delete_int_row(x, i);
    auto fit = ms::estimate_counts(
        reduced, *weights, ms::CountWeighting::unit_weighted);
    REQUIRE(fit.has_value());
    CHECK((summary->coefficient.leave_one_out.row(i).transpose() - fit->estimates)
              .cwiseAbs().maxCoeff() < tol);
  }
}

TEST_CASE("generic jackknife covariance uses the delete-one scaling") {
  ms::RealMat leave_one_out(4, 2);
  leave_one_out << 1.0, 2.0,
                   2.0, 1.0,
                   3.0, 4.0,
                   4.0, 3.0;
  auto covariance = ms::jackknife_covariance(leave_one_out);
  REQUIRE(covariance.has_value());

  const ms::RealVec mean = leave_one_out.colwise().mean().transpose();
  const ms::RealMat centered = leave_one_out.rowwise() - mean.transpose();
  const ms::RealMat expected = 0.75 * centered.transpose() * centered;
  CHECK((*covariance - expected).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("delete-one summaries reject missing data and too-small samples") {
  ms::IntMat x(3, 2);
  x << 0, 0,
       1, ms::na_code,
       0, 1;
  auto weights = ms::loss::identity_weights(2);
  REQUIRE(weights.has_value());
  auto missing = ms::delete_one_kappa_complete(x, *weights);
  REQUIRE(!missing.has_value());
  CHECK(missing.error() == ms::Error::invalid_argument);

  ms::RealMat y(2, 2);
  y << 1.0, 1.0,
       2.0, 2.0;
  auto loss = ms::loss::quadratic_loss(1.0, 2.0);
  REQUIRE(loss.has_value());
  auto small = ms::delete_one_kappa_complete(y, *loss);
  REQUIRE(!small.has_value());
  CHECK(small.error() == ms::Error::invalid_argument);
}
