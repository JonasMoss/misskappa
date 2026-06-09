#include "doctest.h"

#include "misskappa/diagnostics.hpp"
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

TEST_CASE("diagnose_fiml_louis: reports a Louis spectrum on the MAR fixture") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  REQUIRE(W.has_value());

  auto d = ms::diagnose_fiml_louis(x, *W, EmOptions{});
  REQUIRE(d.has_value());

  CHECK(d->c == 3);
  CHECK(d->R == 3);
  CHECK(d->n_subjects > 0u);
  CHECK(d->n_patterns > 0u);

  // Conger's kappa from the diagnostic matches the FIML estimator's Conger.
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());
  CHECK(std::abs(d->kappa_conger - r->estimates(0)) < 1e-6);

  // Spectrum shape: eigenvalues descending; the parallel arrays line up.
  REQUIRE(d->eigenvalues.size() == d->gradient_projection.size());
  REQUIRE(d->eigenvalues.size() == d->variance_contribution.size());
  REQUIRE(d->eigenvalues.size() > 0);
  CHECK(std::abs(d->lambda_max - d->eigenvalues(0)) < tol);
  for (Eigen::Index i = 1; i < d->eigenvalues.size(); ++i) {
    CHECK(d->eigenvalues(i) <= d->eigenvalues(i - 1) + tol);
  }
  CHECK(d->variance >= -tol);
  CHECK(d->retained_rank >= 0);
  CHECK(d->retained_rank <= static_cast<int>(d->eigenvalues.size()));

  // A positive rcond prunes small eigenvalues: positive threshold, no more
  // retained directions than the unpruned run.
  EmOptions pruned{};
  pruned.info_rcond = 0.5;
  auto d_pruned = ms::diagnose_fiml_louis(x, *W, pruned);
  REQUIRE(d_pruned.has_value());
  CHECK(d_pruned->threshold > 0.0);
  CHECK(d_pruned->retained_rank <= d->retained_rank);

  // rcond of zero disables pruning: threshold is exactly zero.
  EmOptions keep_all{};
  keep_all.info_rcond = 0.0;
  auto d_keep = ms::diagnose_fiml_louis(x, *W, keep_all);
  REQUIRE(d_keep.has_value());
  CHECK(d_keep->threshold == 0.0);
}

TEST_CASE("diagnose_fiml_grouped_jackknife: full fit and correction algebra") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  REQUIRE(W.has_value());

  EmOptions opts{};
  opts.tol = 1e-9;
  auto hot = ms::diagnose_fiml_grouped_jackknife(x, *W, opts, 3, true);
  auto cold = ms::diagnose_fiml_grouped_jackknife(x, *W, opts, 3, false);
  REQUIRE(hot.has_value());
  REQUIRE(cold.has_value());

  auto full = ms::estimate_fiml(x, *W, opts);
  REQUIRE(full.has_value());
  CHECK((hot->full_estimates - full->estimates).cwiseAbs().maxCoeff() < 1e-8);
  CHECK((hot->full_vcov - full->vcov).cwiseAbs().maxCoeff() < 1e-8);

  REQUIRE(hot->delete_estimates.rows() == 3);
  REQUIRE(hot->delete_estimates.cols() == 3);
  REQUIRE(hot->delete_iterations.size() == 3);
  CHECK((hot->full_estimates - cold->full_estimates).cwiseAbs().maxCoeff() < 1e-12);
  CHECK((hot->full_vcov - cold->full_vcov).cwiseAbs().maxCoeff() < 1e-12);

  const RealVec delete_mean = cold->delete_estimates.colwise().mean().transpose();
  const RealVec expected_bias = 2.0 * (delete_mean - cold->full_estimates);
  CHECK((cold->jackknife_bias - expected_bias).cwiseAbs().maxCoeff() < 1e-12);
  CHECK((cold->corrected_estimates - (cold->full_estimates - cold->jackknife_bias))
            .cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("diagnose_fiml_penalized: zero penalty and grouped variance") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  REQUIRE(W.has_value());

  EmOptions opts{};
  opts.tol = 1e-9;
  auto full = ms::estimate_fiml(x, *W, opts);
  auto zero = ms::diagnose_fiml_penalized(
      x, *W, opts, ms::FimlPenaltyTarget::uniform, 0.0, 0);
  REQUIRE(full.has_value());
  REQUIRE(zero.has_value());
  CHECK((zero->estimates - full->estimates).cwiseAbs().maxCoeff() < 1e-8);
  CHECK(zero->vcov.rows() == 3);
  CHECK(zero->vcov.cols() == 3);
  CHECK(zero->refits == 0);

  auto penalized = ms::diagnose_fiml_penalized(
      x, *W, opts, ms::FimlPenaltyTarget::independence, 2.0, 3);
  REQUIRE(penalized.has_value());
  CHECK(penalized->groups == 3);
  CHECK(penalized->refits == 3);
  CHECK(penalized->delete_estimates.rows() == 3);
  CHECK(penalized->delete_estimates.cols() == 3);
  CHECK(penalized->delete_iterations.size() == 3);
  CHECK(penalized->estimates.array().isFinite().all());
  CHECK((penalized->vcov - penalized->vcov.transpose()).cwiseAbs().maxCoeff() < 1e-12);
  CHECK(penalized->vcov.diagonal().minCoeff() >= -1e-12);
}

TEST_CASE("diagnose_fiml_louis: rejects degenerate input") {
  auto W = ms::loss::identity_weights(2);
  REQUIRE(W.has_value());

  IntMat one_rater(3, 1);
  one_rater << 0, 1, 0;
  CHECK(ms::diagnose_fiml_louis(one_rater, *W, EmOptions{}).error()
        == ms::Error::invalid_argument);  // R < 2

  RealMat nonsquare(2, 3);
  nonsquare.setZero();
  CHECK(ms::diagnose_fiml_louis(ten_subject_2rater(), nonsquare, EmOptions{}).error()
        == ms::Error::dimension_mismatch);  // weights not square

  IntMat out_of_range(3, 2);
  out_of_range << 0, 1,
                  5, 0,
                  0, 1;
  CHECK(!ms::diagnose_fiml_louis(out_of_range, *W, EmOptions{}).has_value());
}
