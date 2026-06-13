#include "doctest.h"

#include "misskappa/diagnostics.hpp"
#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <vector>

using misskappa::EmOptions;
using misskappa::IntMat;
using misskappa::na_code;
using misskappa::RealMat;
using misskappa::RealVec;
namespace ms = misskappa;

namespace {

constexpr double tol = 1e-9;

void check_three_coef_variance_contract(const ms::Estimation& e, int n) {
  REQUIRE(e.vcov.rows() == 3);
  REQUIRE(e.vcov.cols() == 3);
  CHECK(e.vcov.array().isFinite().all());
  CHECK((e.vcov - e.vcov.transpose()).cwiseAbs().maxCoeff() < 1e-10);
  Eigen::SelfAdjointEigenSolver<RealMat> es(e.vcov);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-8);

  REQUIRE(e.psi.rows() == n);
  REQUIRE(e.psi.cols() == 3);
  CHECK(e.psi.array().isFinite().all());
  const RealMat psi_vcov =
      (e.psi.transpose() * e.psi) / std::pow(static_cast<double>(n), 2);
  CHECK((psi_vcov - e.vcov).cwiseAbs().maxCoeff() < 1e-10);

  REQUIRE(e.null_frac.size() == 3);
  CHECK(e.null_frac.array().isFinite().all());
  CHECK(e.null_frac.minCoeff() >= 0.0);
  CHECK(e.null_frac.maxCoeff() <= 1.0);
}

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

IntMat complete_3rater_2cat() {
  IntMat x(12, 3);
  x <<
    0, 0, 0,
    0, 0, 0,
    0, 0, 1,
    0, 1, 1,
    1, 0, 1,
    1, 1, 0,
    1, 1, 1,
    1, 1, 1,
    0, 1, 0,
    1, 0, 0,
    0, 1, 1,
    1, 0, 1;
  return x;
}

IntMat identified_missing_3rater_3cat() {
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

IntMat frechet_fixture_int() {
  IntMat x(4, 5);
  x << 0, 0, 1, 0, 0,
       0, 1, 2, 1, 1,
       1, 0, 0, 0, 0,
       1, 2, 3, 3, 4;
  return x;
}

IntMat pairwise_only_3rater_2cat() {
  IntMat x(36, 3);
  Eigen::Index row = 0;
  for (int rep = 0; rep < 3; ++rep) {
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) {
        x(row, 0) = a;
        x(row, 1) = b;
        x(row, 2) = na_code;
        ++row;
        x(row, 0) = a;
        x(row, 1) = na_code;
        x(row, 2) = b;
        ++row;
        x(row, 0) = na_code;
        x(row, 1) = a;
        x(row, 2) = b;
        ++row;
      }
    }
  }
  return x;
}

IntMat chain_missing_pair_3rater_2cat() {
  IntMat x(12, 3);
  x << 0, 0, na_code,
       0, 0, na_code,
       0, 0, na_code,
       1, 1, na_code,
       1, 1, na_code,
       0, 1, na_code,
       na_code, 0, 0,
       na_code, 0, 0,
       na_code, 1, 1,
       na_code, 1, 1,
       na_code, 1, 0,
       na_code, 0, 1;
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

TEST_CASE("estimate_fiml: complete data supports non-symmetric weights") {
  IntMat x = complete_3rater_2cat();
  RealMat W(2, 2);
  W << 1.0, 0.2,
       0.7, 1.0;

  auto available = ms::estimate_available(x, W);
  auto fiml = ms::estimate_fiml(x, W, EmOptions{});
  REQUIRE(available.has_value());
  REQUIRE(fiml.has_value());

  CHECK((fiml->estimates - available->estimates).cwiseAbs().maxCoeff() < 1e-10);
  check_three_coef_variance_contract(*fiml, static_cast<int>(x.rows()));
}

TEST_CASE("estimate_fiml: sparse MAR fixture succeeds with null-frac diagnostic") {
  // All rater pairs are co-observed, so the coefficients are estimable even
  // though the saturated 3^3 nuisance is far from identified at n = 12. The
  // former hard not_identified gate is now the null_frac diagnostic.
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());

  CHECK(r->estimates.array().isFinite().all());
  CHECK(r->vcov.array().isFinite().all());
  REQUIRE(r->null_frac.size() == 3);
  CHECK(r->null_frac.minCoeff() >= 0.0);
  CHECK(r->null_frac.maxCoeff() <= 1.0);
  // Sample information is rank-deficient along directions the coefficients
  // touch; the diagnostic must say so.
  CHECK(r->null_frac.maxCoeff() > 1e-3);
}

TEST_CASE("estimate_fiml: complete data has zero null-frac diagnostic") {
  IntMat x = ten_subject_2rater();
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());
  REQUIRE(r->null_frac.size() == 3);
  CHECK(r->null_frac.maxCoeff() < 1e-6);
}

TEST_CASE("estimate_fiml: flatten selects a unique nearby posterior mode") {
  IntMat x = twelve_subject_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);

  auto strict = ms::estimate_fiml(x, *W, EmOptions{});
  EmOptions flat_opts;
  flat_opts.flatten = 0.1;
  auto flat = ms::estimate_fiml(x, *W, flat_opts);
  REQUIRE(strict.has_value());
  REQUIRE(flat.has_value());

  CHECK(flat->estimates.array().isFinite().all());
  CHECK(flat->vcov.array().isFinite().all());
  // Total pseudo-mass 0.1 against n = 12 shrinks the fitted table toward
  // uniform with weight ~ 0.1 / 12.1, so the kappas move a little but stay
  // close to the strict-ML face.
  CHECK((flat->estimates - strict->estimates).cwiseAbs().maxCoeff() < 0.1);

  // The flattened start is independent of start_alpha: the posterior mode is
  // unique, so two different starts agree far beyond face width.
  EmOptions other_start = flat_opts;
  other_start.start_alpha = 1.0;
  auto flat2 = ms::estimate_fiml(x, *W, other_start);
  REQUIRE(flat2.has_value());
  CHECK((flat->estimates - flat2->estimates).cwiseAbs().maxCoeff() < 1e-5);
}

TEST_CASE("estimate_fiml: identified missing fixture returns finite estimates") {
  IntMat x = identified_missing_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());

  CHECK(r->estimates.array().isFinite().all());
  CHECK(r->vcov.array().isFinite().all());
  CHECK(r->vcov.diagonal().minCoeff() >= -1e-10);
}

TEST_CASE("estimate_fiml: variance is symmetric and PSD") {
  IntMat x = identified_missing_3rater_3cat();
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
  IntMat x = identified_missing_3rater_3cat();
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

TEST_CASE("estimate_fiml: benign nuisance non-identification succeeds") {
  IntMat x = pairwise_only_3rater_2cat();
  auto W = ms::loss::identity_weights(2);
  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(r.has_value());

  CHECK(r->estimates.size() == 3);
  CHECK(r->estimates.array().isFinite().all());
  CHECK(r->vcov.array().isFinite().all());
  CHECK((r->vcov - r->vcov.transpose()).cwiseAbs().maxCoeff() < 1e-10);

  Eigen::SelfAdjointEigenSolver<RealMat> es(r->vcov);
  REQUIRE(es.info() == Eigen::Success);
  CHECK(es.eigenvalues().minCoeff() > -1e-8);

  REQUIRE(r->psi.rows() == x.rows());
  REQUIRE(r->psi.cols() == 3);
  const RealMat psi_vcov =
      (r->psi.transpose() * r->psi) / std::pow(static_cast<double>(x.rows()), 2);
  CHECK((psi_vcov - r->vcov).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE("estimate_fiml: missing rater-pair coefficient non-identification errors") {
  IntMat x = chain_missing_pair_3rater_2cat();
  auto W = ms::loss::identity_weights(2);

  auto r = ms::estimate_fiml(x, *W, EmOptions{});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::not_identified);

  const std::vector<RealMat> weights{*W};
  auto many = ms::estimate_fiml_many(x, weights, EmOptions{});
  REQUIRE(!many.has_value());
  CHECK(many.error() == ms::Error::not_identified);
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
  IntMat x = identified_missing_3rater_3cat();
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

TEST_CASE("estimate_fiml_gwise: missing requested rater tuple errors") {
  IntMat x = pairwise_only_3rater_2cat();
  auto distance = ms::loss::frechet_nominal_distance(2);
  REQUIRE(distance.has_value());

  auto r = ms::estimate_fiml_gwise(x, *distance, EmOptions{}, ms::GwiseOptions{3});
  REQUIRE(!r.has_value());
  CHECK(r.error() == ms::Error::not_identified);
}

TEST_CASE("estimate_fiml_gwise: variance is symmetric PSD and psi reconstructs it") {
  IntMat x = identified_missing_3rater_3cat();
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
  IntMat x = identified_missing_3rater_3cat();
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

TEST_CASE("estimate_fiml: C++ SQUAREM acceleration preserves estimates") {
  IntMat x = identified_missing_3rater_3cat();
  RealVec v(3);
  v << 1.0, 2.0, 3.0;
  auto W = ms::loss::quadratic_weights(3, v);

  EmOptions plain{};
  plain.tol = 1e-9;
  plain.max_iter = 20000;
  auto r_plain = ms::estimate_fiml(x, *W, plain);
  REQUIRE(r_plain.has_value());

  EmOptions squarem = plain;
  squarem.acceleration = ms::EmAcceleration::squarem;
  auto r_squarem = ms::estimate_fiml(x, *W, squarem);
  REQUIRE(r_squarem.has_value());

  CHECK((r_plain->estimates - r_squarem->estimates).cwiseAbs().maxCoeff() < 1e-6);
  CHECK(r_squarem->vcov.allFinite());
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
  IntMat x = identified_missing_3rater_3cat();
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
