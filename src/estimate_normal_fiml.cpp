#include "detail_pattern_checks.hpp"
#include "misskappa/estimate.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace misskappa {

namespace {

constexpr double inverse_rcond_tol = 1e-12;
constexpr double zero_tol = 1e-9;

struct Pattern {
  std::vector<int> rows;
  std::vector<int> observed;
  std::vector<int> missing;
  RealMat x_observed;
  RealVec sum_observed;
  RealMat cross_observed;
};

struct NormalEmFit {
  RealVec mu;
  RealMat sigma;
  std::vector<Pattern> patterns;
  RealMat x;
  int iterations = 0;
  bool converged = false;
};

struct GradientResult {
  RealVec estimates;
  RealMat gradient;  // K x q, q = p + vech_size(p).
};

int vech_size(int p) {
  return p * (p + 1) / 2;
}

RealMat compact_observed_rows(RealMatView ratings) {
  int keep = 0;
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    bool any = false;
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      if (std::isfinite(ratings(i, j))) {
        any = true;
        break;
      }
    }
    if (any) ++keep;
  }

  RealMat out(keep, ratings.cols());
  int row = 0;
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    bool any = false;
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      if (std::isfinite(ratings(i, j))) {
        any = true;
        break;
      }
    }
    if (!any) continue;
    out.row(row) = ratings.row(i);
    ++row;
  }
  return out;
}

Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> vech_positions(int p) {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> pos(p, p);
  pos.setZero();
  int k = 0;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      pos(i, j) = k;
      pos(j, i) = k;
      ++k;
    }
  }
  return pos;
}

RealVec gather_vector(const RealVec& x, const std::vector<int>& idx) {
  RealVec out(static_cast<Eigen::Index>(idx.size()));
  for (std::size_t k = 0; k < idx.size(); ++k) out(static_cast<Eigen::Index>(k)) = x(idx[k]);
  return out;
}

RealMat gather_matrix(
    const RealMat& x, const std::vector<int>& rows, const std::vector<int>& cols) {
  RealMat out(static_cast<Eigen::Index>(rows.size()), static_cast<Eigen::Index>(cols.size()));
  for (std::size_t i = 0; i < rows.size(); ++i) {
    for (std::size_t j = 0; j < cols.size(); ++j) {
      out(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = x(rows[i], cols[j]);
    }
  }
  return out;
}

RealMat gather_square(const RealMat& x, const std::vector<int>& idx) {
  return gather_matrix(x, idx, idx);
}

RealMat gather_submatrix(
    const RealMat& x, const std::vector<int>& rows, const std::vector<int>& cols) {
  return gather_matrix(x, rows, cols);
}

Result<RealMat> inverse_matrix(const RealMat& x) {
  if (x.rows() != x.cols() || x.rows() == 0) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  if (!x.allFinite()) return misskappa::unexpected(Error::numerical_error);

  Eigen::LLT<RealMat> llt(x);
  if (llt.info() == Eigen::Success) {
    const double rcond = llt.rcond();
    if (!std::isfinite(rcond) || rcond <= inverse_rcond_tol) {
      return misskappa::unexpected(Error::numerical_error);
    }
    RealMat inv = llt.solve(RealMat::Identity(x.rows(), x.cols()));
    if (!inv.allFinite()) return misskappa::unexpected(Error::numerical_error);
    inv = 0.5 * (inv + inv.transpose());
    return inv;
  }

  Eigen::FullPivLU<RealMat> lu(x);
  lu.setThreshold(inverse_rcond_tol);
  if (!lu.isInvertible()) return misskappa::unexpected(Error::numerical_error);
  RealMat inv = lu.inverse();
  if (!inv.allFinite()) return misskappa::unexpected(Error::numerical_error);
  return inv;
}

Result<RealMat> solve_matrix(const RealMat& a, const RealMat& b) {
  Eigen::FullPivLU<RealMat> lu(a);
  lu.setThreshold(1e-12);
  if (!lu.isInvertible()) return misskappa::unexpected(Error::numerical_error);
  RealMat x = lu.solve(b);
  if (!x.allFinite()) return misskappa::unexpected(Error::numerical_error);
  return x;
}

Result<double> quadratic_bp_chance(RealMatView ratings, const RealVec& values) {
  if (values.size() > 0) {
    if (!values.allFinite()) return misskappa::unexpected(Error::invalid_argument);
    const double C = static_cast<double>(values.size());
    const double chance =
        (2.0 / (C * C)) * (C * values.squaredNorm() - std::pow(values.sum(), 2));
    if (!std::isfinite(chance) || chance <= zero_tol) {
      return misskappa::unexpected(Error::invalid_argument);
    }
    return chance;
  }

  bool any = false;
  double lo = std::numeric_limits<double>::infinity();
  double hi = -std::numeric_limits<double>::infinity();
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      const double x = ratings(i, j);
      if (!std::isfinite(x)) continue;
      any = true;
      lo = std::min(lo, x);
      hi = std::max(hi, x);
    }
  }
  if (!any || !(hi > lo)) return misskappa::unexpected(Error::invalid_argument);
  const double chance = (hi - lo) * (hi - lo) / 6.0;
  if (!std::isfinite(chance) || chance <= zero_tol) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  return chance;
}

// Project a p x p saturated covariance onto the scalar compound symmetry
// Sigma = psi J + e I (the block_size == 1 case).
Result<RealMat> project_scalar_variance_components(const RealMat& sigma) {
  const int p = static_cast<int>(sigma.rows());
  if (p < 2 || sigma.cols() != p || !sigma.allFinite()) {
    return misskappa::unexpected(Error::numerical_error);
  }
  const double trace = sigma.diagonal().sum();
  const double total = sigma.sum();
  const double contrast_eigen =
      (trace - total / static_cast<double>(p)) / static_cast<double>(p - 1);
  const double mean_eigen = total / static_cast<double>(p);

  double e = 0.0;
  double psi = 0.0;
  if (mean_eigen >= contrast_eigen) {
    e = contrast_eigen;
    psi = (mean_eigen - contrast_eigen) / static_cast<double>(p);
  } else {
    e = trace / static_cast<double>(p);
    psi = 0.0;
  }

  const double lambda_mean = e + static_cast<double>(p) * psi;
  constexpr double min_eigen = 1e-12;
  if (!std::isfinite(e) || !std::isfinite(psi) || !std::isfinite(lambda_mean) ||
      e <= min_eigen || lambda_mean <= min_eigen || psi < -1e-12) {
    return misskappa::unexpected(Error::numerical_error);
  }

  RealMat out = RealMat::Constant(p, p, std::max(0.0, psi));
  out.diagonal().array() += e;
  return out;
}

// Project a symmetric d x d matrix onto the PSD cone with a floor on the
// eigenvalues. Returns the unchanged matrix when it is already interior.
RealMat clip_psd(const RealMat& M, double floor) {
  Eigen::SelfAdjointEigenSolver<RealMat> es(0.5 * (M + M.transpose()));
  RealVec lam = es.eigenvalues();
  if (lam.minCoeff() >= floor) return M;
  for (int i = 0; i < lam.size(); ++i) lam(i) = std::max(lam(i), floor);
  return es.eigenvectors() * lam.asDiagonal() * es.eigenvectors().transpose();
}

// Project a (R*d) x (R*d) saturated covariance onto the block compound symmetry
// Sigma = J_R (x) Psi + I_R (x) E. The constrained Gaussian MLE is block
// averaging (the rater eigenbasis simultaneously diagonalises every candidate
// Sigma, so the likelihood factorises): Psi = mean off-diagonal d x d block,
// E = mean diagonal block - Psi. PSD is guarded through the two eigen-blocks E
// (contrast directions) and R*Psi + E (mean direction); interior fits are left
// exactly at the block average.
Result<RealMat> project_block_variance_components(const RealMat& sigma, int d) {
  const int p = static_cast<int>(sigma.rows());
  if (d < 1 || p < 2 * d || p % d != 0 || sigma.cols() != p || !sigma.allFinite()) {
    return misskappa::unexpected(Error::numerical_error);
  }
  const int R = p / d;
  RealMat diag_sum = RealMat::Zero(d, d);
  RealMat off_sum = RealMat::Zero(d, d);
  for (int r = 0; r < R; ++r) {
    diag_sum += sigma.block(r * d, r * d, d, d);
    for (int s = 0; s < R; ++s) {
      if (s == r) continue;
      off_sum += sigma.block(r * d, s * d, d, d);
    }
  }
  RealMat Dbar = diag_sum / static_cast<double>(R);
  RealMat Psi = off_sum / static_cast<double>(R * (R - 1));
  Psi = 0.5 * (Psi + Psi.transpose());
  RealMat E = Dbar - Psi;
  E = 0.5 * (E + E.transpose());

  constexpr double min_eigen = 1e-12;
  RealMat Mmean = static_cast<double>(R) * Psi + E;          // mean-direction block
  Eigen::SelfAdjointEigenSolver<RealMat> esE(E);
  Eigen::SelfAdjointEigenSolver<RealMat> esM(Mmean);
  if (esE.info() != Eigen::Success || esM.info() != Eigen::Success) {
    return misskappa::unexpected(Error::numerical_error);
  }
  if (esE.eigenvalues().minCoeff() < min_eigen ||
      esM.eigenvalues().minCoeff() < min_eigen) {
    E = clip_psd(E, min_eigen);
    Mmean = clip_psd(Mmean, min_eigen);
    Psi = (Mmean - E) / static_cast<double>(R);
    Psi = 0.5 * (Psi + Psi.transpose());
  }

  RealMat out(p, p);
  for (int r = 0; r < R; ++r) {
    for (int s = 0; s < R; ++s) {
      out.block(r * d, s * d, d, d) = (r == s) ? (Psi + E).eval() : Psi;
    }
  }
  if (!out.allFinite()) return misskappa::unexpected(Error::numerical_error);
  return out;
}

Result<RealMat> project_variance_components(const RealMat& sigma, int block_size = 1) {
  if (block_size <= 1) return project_scalar_variance_components(sigma);
  return project_block_variance_components(sigma, block_size);
}

// Jacobian d(mu, vech Sigma) / d theta for the variance-components model, where
// theta = (mu params, vech Psi, vech E). block_size == 1 is the scalar case
// (Psi, E scalar; one psi and one e column). For block_size = d > 1 the columns
// are the d x d blocks Psi and E in vech order, and Sigma is J_R (x) Psi +
// I_R (x) E with rater-major stacked columns.
RealMat variance_component_jacobian(int p, NormalMeanModel mean_model,
                                    int block_size = 1) {
  const int d = (block_size <= 1) ? 1 : block_size;
  const int R = p / d;
  const int sd = vech_size(d);                 // params per d x d block
  const int q = p + vech_size(p);
  const int mean_params = (mean_model == NormalMeanModel::equal) ? d : p;
  RealMat A = RealMat::Zero(q, mean_params + 2 * sd);

  if (mean_model == NormalMeanModel::equal) {
    // shared d-vector mean repeated across raters
    for (int r = 0; r < R; ++r) {
      for (int a = 0; a < d; ++a) A(r * d + a, a) = 1.0;
    }
  } else {
    A.block(0, 0, p, p).setIdentity();
  }

  // d x d vech index for a feature pair (column-major lower triangle), matching
  // the (j outer, i >= j inner) order used for the full vech below.
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> vd(d, d);
  {
    int t = 0;
    for (int b = 0; b < d; ++b) {
      for (int a = b; a < d; ++a) { vd(a, b) = t; vd(b, a) = t; ++t; }
    }
  }

  const int psi0 = mean_params;
  const int e0 = mean_params + sd;
  int k = 0;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      const int fa = i % d, fb = j % d;
      A(p + k, psi0 + vd(fa, fb)) = 1.0;                 // d Sigma_ij / d Psi
      if (i / d == j / d) A(p + k, e0 + vd(fa, fb)) = 1.0;  // same rater -> + E
      ++k;
    }
  }
  return A;
}

Result<std::vector<Pattern>> build_patterns(const RealMat& x) {
  const int n = static_cast<int>(x.rows());
  const int p = static_cast<int>(x.cols());
  std::map<std::string, std::vector<int>> groups;
  for (int i = 0; i < n; ++i) {
    std::string key;
    key.reserve(static_cast<std::size_t>(p));
    bool any = false;
    for (int j = 0; j < p; ++j) {
      const bool obs = std::isfinite(x(i, j));
      key.push_back(obs ? '1' : '0');
      any = any || obs;
    }
    if (!any) return misskappa::unexpected(Error::invalid_argument);
    groups[key].push_back(i);
  }

  std::vector<Pattern> patterns;
  patterns.reserve(groups.size());
  for (const auto& kv : groups) {
    Pattern pat;
    pat.rows = kv.second;
    for (int j = 0; j < p; ++j) {
      if (kv.first[static_cast<std::size_t>(j)] == '1') {
        pat.observed.push_back(j);
      } else {
        pat.missing.push_back(j);
      }
    }
    if (pat.observed.empty()) return misskappa::unexpected(Error::invalid_argument);
    pat.x_observed = gather_matrix(x, pat.rows, pat.observed);
    pat.sum_observed = pat.x_observed.colwise().sum().transpose();
    pat.cross_observed = pat.x_observed.transpose() * pat.x_observed;
    patterns.push_back(std::move(pat));
  }
  return patterns;
}

Result<NormalEmFit> fit_saturated_normal(RealMatView ratings, EmOptions opts) {
  if (ratings.rows() < 1 || ratings.cols() < 2) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  if (opts.max_iter < 1 || !(opts.tol > 0.0) || !std::isfinite(opts.tol)) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  if (opts.mean == NormalMeanModel::equal &&
      opts.covariance != NormalCovarianceModel::variance_components) {
    return misskappa::unexpected(Error::invalid_argument);
  }

  NormalEmFit out;
  out.x = compact_observed_rows(ratings);
  const int n = static_cast<int>(out.x.rows());
  const int p = static_cast<int>(out.x.cols());
  if (n < 1 || p < 2) return misskappa::unexpected(Error::invalid_argument);

  const auto mask = detail::finite_mask(out.x);
  auto identified = (opts.covariance == NormalCovarianceModel::variance_components)
      ? detail::require_connected_observation(mask)
      : detail::require_complete_pair_observation(mask);
  if (!identified) return misskappa::unexpected(identified.error());

  auto patterns = build_patterns(out.x);
  if (!patterns) return misskappa::unexpected(patterns.error());
  out.patterns = std::move(*patterns);

  RealVec count = RealVec::Zero(p);
  RealVec sum = RealVec::Zero(p);
  RealVec sumsq = RealVec::Zero(p);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      const double v = out.x(i, j);
      if (!std::isfinite(v)) continue;
      count(j) += 1.0;
      sum(j) += v;
      sumsq(j) += v * v;
    }
  }

  out.mu.resize(p);
  RealVec variance(p);
  for (int j = 0; j < p; ++j) {
    if (count(j) <= 0.0) return misskappa::unexpected(Error::not_identified);
    out.mu(j) = sum(j) / count(j);
    double v = 1.0;
    if (count(j) > 1.0) {
      v = (sumsq(j) - sum(j) * sum(j) / count(j)) / (count(j) - 1.0);
    }
    if (!std::isfinite(v) || v <= 0.0) v = 1.0;
    variance(j) = v;
  }
  if (opts.mean == NormalMeanModel::equal) {
    out.mu.setConstant(sum.sum() / count.sum());
  }
  out.sigma = variance.asDiagonal();

  for (int iter = 1; iter <= opts.max_iter; ++iter) {
    RealVec t1 = RealVec::Zero(p);
    RealMat t2 = RealMat::Zero(p, p);

    for (const Pattern& pat : out.patterns) {
      const double ng = static_cast<double>(pat.rows.size());
      const int no = static_cast<int>(pat.observed.size());
      for (int a = 0; a < no; ++a) {
        const int ja = pat.observed[static_cast<std::size_t>(a)];
        t1(ja) += pat.sum_observed(a);
        for (int b = 0; b < no; ++b) {
          const int jb = pat.observed[static_cast<std::size_t>(b)];
          t2(ja, jb) += pat.cross_observed(a, b);
        }
      }

      if (!pat.missing.empty()) {
        auto soo_inv = inverse_matrix(gather_square(out.sigma, pat.observed));
        if (!soo_inv) return misskappa::unexpected(soo_inv.error());
        const RealMat smo = gather_submatrix(out.sigma, pat.missing, pat.observed);
        const RealMat smm = gather_square(out.sigma, pat.missing);
        const RealMat b = smo * (*soo_inv);
        const RealVec mu_o = gather_vector(out.mu, pat.observed);
        const RealVec mu_m = gather_vector(out.mu, pat.missing);

        const RealVec se = pat.sum_observed - ng * mu_o;
        const RealMat see = pat.cross_observed - pat.sum_observed * mu_o.transpose() -
                            mu_o * pat.sum_observed.transpose() + ng * mu_o * mu_o.transpose();

        const RealVec sum_missing = b * se + ng * mu_m;
        const RealMat xo_e_cross = pat.cross_observed - pat.sum_observed * mu_o.transpose();
        const RealMat cross_observed_missing =
            xo_e_cross * b.transpose() + pat.sum_observed * mu_m.transpose();
        const RealMat b_see_b = b * see * b.transpose();
        const RealVec b_se = b * se;
        const RealMat cross_missing =
            b_see_b + b_se * mu_m.transpose() + mu_m * b_se.transpose() +
            ng * mu_m * mu_m.transpose();
        const RealMat conditional_cov = smm - b * smo.transpose();

        const int nm = static_cast<int>(pat.missing.size());
        for (int a = 0; a < nm; ++a) {
          const int ja = pat.missing[static_cast<std::size_t>(a)];
          t1(ja) += sum_missing(a);
          for (int c = 0; c < no; ++c) {
            const int jc = pat.observed[static_cast<std::size_t>(c)];
            const double value = cross_observed_missing(c, a);
            t2(jc, ja) += value;
            t2(ja, jc) += value;
          }
          for (int c = 0; c < nm; ++c) {
            const int jc = pat.missing[static_cast<std::size_t>(c)];
            t2(ja, jc) += cross_missing(a, c) + ng * conditional_cov(a, c);
          }
        }
      }
    }

    RealVec mu_new = t1 / static_cast<double>(n);
    if (opts.mean == NormalMeanModel::equal) {
      mu_new.setConstant(t1.sum() / (static_cast<double>(n) * static_cast<double>(p)));
    }
    RealMat sigma_new = t2 / static_cast<double>(n) - mu_new * mu_new.transpose();
    sigma_new = 0.5 * (sigma_new + sigma_new.transpose());
    if (opts.covariance == NormalCovarianceModel::variance_components) {
      auto projected = project_variance_components(sigma_new, opts.block_size);
      if (!projected) return misskappa::unexpected(projected.error());
      sigma_new = std::move(*projected);
    }
    const double delta = std::max(
        (mu_new - out.mu).cwiseAbs().maxCoeff(),
        (sigma_new - out.sigma).cwiseAbs().maxCoeff());
    out.mu = std::move(mu_new);
    out.sigma = std::move(sigma_new);
    out.iterations = iter;
    if (delta < opts.tol) {
      out.converged = true;
      break;
    }
  }

  if (!out.mu.allFinite() || !out.sigma.allFinite()) {
    return misskappa::unexpected(Error::numerical_error);
  }
  return out;
}

Result<RealMat> score_matrix(
    const RealVec& mu, const RealMat& sigma, const RealMat& x,
    const std::vector<Pattern>& patterns,
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& vech_pos) {
  const int n = static_cast<int>(x.rows());
  const int p = static_cast<int>(mu.size());
  const int q = p + vech_size(p);
  RealMat scores = RealMat::Zero(n, q);

  for (const Pattern& pat : patterns) {
    const int no = static_cast<int>(pat.observed.size());
    auto soo_inv = inverse_matrix(gather_square(sigma, pat.observed));
    if (!soo_inv) return misskappa::unexpected(soo_inv.error());
    const RealVec mu_o = gather_vector(mu, pat.observed);
    RealMat e = pat.x_observed.rowwise() - mu_o.transpose();
    RealMat u = e * (*soo_inv);
    const RealMat half_const = 0.5 * (*soo_inv);

    for (std::size_t rr = 0; rr < pat.rows.size(); ++rr) {
      const int row = pat.rows[rr];
      for (int a = 0; a < no; ++a) scores(row, pat.observed[static_cast<std::size_t>(a)]) += u(rr, a);
      for (int a = 0; a < no; ++a) {
        for (int b = 0; b <= a; ++b) {
          const int ja = pat.observed[static_cast<std::size_t>(a)];
          const int jb = pat.observed[static_cast<std::size_t>(b)];
          const double wt = (ja == jb) ? 1.0 : 2.0;
          const int col = p + vech_pos(ja, jb);
          scores(row, col) += wt * (0.5 * u(rr, a) * u(rr, b) - half_const(a, b));
        }
      }
    }
  }

  return scores;
}

Result<RealMat> observed_information(
    const RealVec& mu, const RealMat& sigma, const std::vector<Pattern>& patterns,
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& vech_pos, int n, int p) {
  const int q = p + vech_size(p);
  RealMat H = RealMat::Zero(q, q);

  for (const Pattern& pat : patterns) {
    const int no = static_cast<int>(pat.observed.size());
    const double ng = static_cast<double>(pat.rows.size());
    auto soo_inv = inverse_matrix(gather_square(sigma, pat.observed));
    if (!soo_inv) return misskappa::unexpected(soo_inv.error());
    const RealMat& a_inv = *soo_inv;
    const RealVec mu_o = gather_vector(mu, pat.observed);
    const RealVec se = pat.sum_observed - ng * mu_o;
    const RealMat see = pat.cross_observed - pat.sum_observed * mu_o.transpose() -
                        mu_o * pat.sum_observed.transpose() + ng * mu_o * mu_o.transpose();
    const RealVec v = a_inv * se;
    const RealMat u = a_inv * see * a_inv;

    for (int a = 0; a < no; ++a) {
      const int ja = pat.observed[static_cast<std::size_t>(a)];
      for (int c = 0; c < no; ++c) {
        H(ja, pat.observed[static_cast<std::size_t>(c)]) -= ng * a_inv(a, c);
      }
    }

    for (int a = 0; a < no; ++a) {
      for (int b = 0; b <= a; ++b) {
        const int ja = pat.observed[static_cast<std::size_t>(a)];
        const int jb = pat.observed[static_cast<std::size_t>(b)];
        const int col = p + vech_pos(ja, jb);
        for (int c = 0; c < no; ++c) {
          const int jc = pat.observed[static_cast<std::size_t>(c)];
          const double value = (a == b) ? -a_inv(c, a) * v(a)
                                        : -(a_inv(c, a) * v(b) + a_inv(c, b) * v(a));
          H(jc, col) += value;
          H(col, jc) += value;
        }
      }
    }

    for (int a = 0; a < no; ++a) {
      for (int b = 0; b <= a; ++b) {
        const int ja = pat.observed[static_cast<std::size_t>(a)];
        const int jb = pat.observed[static_cast<std::size_t>(b)];
        const int row = p + vech_pos(ja, jb);
        const double row_weight = (ja == jb) ? 1.0 : 2.0;
        for (int c = 0; c < no; ++c) {
          for (int d = 0; d <= c; ++d) {
            const int jc = pat.observed[static_cast<std::size_t>(c)];
            const int jd = pat.observed[static_cast<std::size_t>(d)];
            const int col = p + vech_pos(jc, jd);

            double aa = 0.0;
            double au = 0.0;
            double ua = 0.0;
            if (c == d) {
              aa = a_inv(a, c) * a_inv(c, b);
              au = a_inv(a, c) * u(c, b);
              ua = u(a, c) * a_inv(c, b);
            } else {
              aa = a_inv(a, c) * a_inv(d, b) + a_inv(a, d) * a_inv(c, b);
              au = a_inv(a, c) * u(d, b) + a_inv(a, d) * u(c, b);
              ua = u(a, c) * a_inv(d, b) + u(a, d) * a_inv(c, b);
            }
            H(row, col) += row_weight * 0.5 * (ng * aa - au - ua);
          }
        }
      }
    }
  }

  H = -H / static_cast<double>(n);
  H = 0.5 * (H + H.transpose());
  if (!H.allFinite()) return misskappa::unexpected(Error::numerical_error);
  return H;
}

template <typename GradientFn>
Result<NormalFimlEstimation> estimate_normal_fiml_impl(
    RealMatView ratings, EmOptions opts, GradientFn gradient_fn) {
  auto em = fit_saturated_normal(ratings, opts);
  if (!em) return misskappa::unexpected(em.error());
  const int n = static_cast<int>(em->x.rows());
  const int p = static_cast<int>(em->x.cols());
  const int pstar = vech_size(p);
  const int q = p + pstar;
  const auto vech_pos = vech_positions(p);

  GradientResult grad = gradient_fn(em->mu, em->sigma);
  if (grad.gradient.cols() != q || grad.gradient.rows() != grad.estimates.size()) {
    return misskappa::unexpected(Error::dimension_mismatch);
  }
  if (!grad.estimates.allFinite() || !grad.gradient.allFinite()) {
    return misskappa::unexpected(Error::numerical_error);
  }

  auto scores = score_matrix(em->mu, em->sigma, em->x, em->patterns, vech_pos);
  if (!scores) return misskappa::unexpected(scores.error());
  auto H = observed_information(em->mu, em->sigma, em->patterns, vech_pos, n, p);
  if (!H) return misskappa::unexpected(H.error());
  Result<RealMat> solved = misskappa::unexpected(Error::numerical_error);
  if (opts.covariance == NormalCovarianceModel::variance_components) {
    const RealMat A = variance_component_jacobian(p, opts.mean, opts.block_size);
    const RealMat H_reduced = A.transpose() * (*H) * A;
    const RealMat G_reduced_t = (grad.gradient * A).transpose();
    auto reduced = solve_matrix(H_reduced, G_reduced_t);
    if (!reduced) return misskappa::unexpected(reduced.error());
    solved = RealMat(A * (*reduced));
  } else {
    solved = solve_matrix(*H, grad.gradient.transpose());
    if (!solved) return misskappa::unexpected(solved.error());
  }

  RealMat psi = (*scores) * (*solved);
  RealMat vcov = (psi.transpose() * psi) / std::pow(static_cast<double>(n), 2);
  vcov = 0.5 * (vcov + vcov.transpose());

  Estimation fit;
  fit.estimates = std::move(grad.estimates);
  fit.vcov = std::move(vcov);
  fit.psi = std::move(psi);
  fit.null_frac = RealVec();

  NormalFimlEstimation out;
  out.fit = std::move(fit);
  out.mu = em->mu;
  out.sigma = em->sigma;
  out.iterations = em->iterations;
  out.converged = em->converged;
  return out;
}

GradientResult alpha_gradient(const RealVec& mu, const RealMat& sigma) {
  const int p = static_cast<int>(mu.size());
  const int pstar = vech_size(p);
  const int q = p + pstar;
  const double trace = sigma.diagonal().sum();
  const double total = sigma.sum();
  const double factor = static_cast<double>(p) / static_cast<double>(p - 1);

  GradientResult out;
  out.estimates = RealVec(1);
  out.estimates(0) = factor * (1.0 - trace / total);
  out.gradient = RealMat::Zero(1, q);

  int k = 0;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      const double a = (i == j) ? 1.0 : 0.0;
      const double b = (i == j) ? 1.0 : 2.0;
      out.gradient(0, p + k) = factor * (trace * b / (total * total) - a / total);
      ++k;
    }
  }
  (void)mu;
  return out;
}

GradientResult scalar_quadratic_gradient(
    const RealVec& mu, const RealMat& sigma, double bp_chance) {
  const int R = static_cast<int>(mu.size());
  const int pstar = vech_size(R);
  const int q = R + pstar;
  const double t1 = sigma.sum();
  const double t2 = sigma.diagonal().sum();
  const double mubar = mu.mean();
  const double t3 = (mu.array() - mubar).square().sum();

  const double dC = static_cast<double>(R - 1) * t2 + static_cast<double>(R) * t3;
  const double nC = t1 - t2;
  const double dF = static_cast<double>(R - 1) * (t2 + t3);
  const double nF = t1 - t2 - t3;
  const double d_obs =
      (2.0 / static_cast<double>(R - 1)) * (t2 + t3 - t1 / static_cast<double>(R));

  RealMat At = RealMat::Zero(3, q);
  int k = 0;
  for (int j = 0; j < R; ++j) {
    for (int i = j; i < R; ++i) {
      const bool diag = (i == j);
      At(0, R + k) = diag ? 1.0 : 2.0;
      At(1, R + k) = diag ? 1.0 : 0.0;
      ++k;
    }
  }
  for (int j = 0; j < R; ++j) At(2, j) = 2.0 * (mu(j) - mubar);

  RealMat Jh = RealMat::Zero(3, 3);
  Jh(0, 0) = 1.0 / dC;
  Jh(0, 1) = -1.0 / dC - nC * static_cast<double>(R - 1) / (dC * dC);
  Jh(0, 2) = -nC * static_cast<double>(R) / (dC * dC);
  const double f_t23 = -1.0 / dF - nF * static_cast<double>(R - 1) / (dF * dF);
  Jh(1, 0) = 1.0 / dF;
  Jh(1, 1) = f_t23;
  Jh(1, 2) = f_t23;
  const double bp_factor = -2.0 / (bp_chance * static_cast<double>(R - 1));
  Jh(2, 0) = -bp_factor / static_cast<double>(R);
  Jh(2, 1) = bp_factor;
  Jh(2, 2) = bp_factor;

  GradientResult out;
  out.estimates = RealVec(3);
  out.estimates(0) = nC / dC;
  out.estimates(1) = nF / dF;
  out.estimates(2) = 1.0 - d_obs / bp_chance;
  out.gradient = Jh * At;
  return out;
}

int vector_col(int rater, int feature, int features) {
  return rater * features + feature;
}

RealMat vector_design(int R, int features, const RealMat& W, int kind) {
  const int q = R * features;
  RealMat out = RealMat::Zero(q, q);
  for (int r = 0; r < R; ++r) {
    for (int s = 0; s < R; ++s) {
      double rater_weight = 1.0;
      if (kind == 1) {
        rater_weight = (r == s) ? 1.0 : 0.0;
      } else if (kind == 2) {
        rater_weight = (r == s ? 1.0 : 0.0) - 1.0 / static_cast<double>(R);
      }
      if (rater_weight == 0.0) continue;
      for (int a = 0; a < features; ++a) {
        for (int b = 0; b < features; ++b) {
          out(vector_col(r, a, features), vector_col(s, b, features)) =
              rater_weight * W(a, b);
        }
      }
    }
  }
  return out;
}

RealVec vech_gradient_from_matrix(const RealMat& A) {
  const int p = static_cast<int>(A.rows());
  RealVec out(vech_size(p));
  int k = 0;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      out(k) = (i == j) ? A(i, j) : A(i, j) + A(j, i);
      ++k;
    }
  }
  return out;
}

Result<void> validate_feature_metric(const RealMat& W) {
  if (W.rows() != W.cols() || W.rows() < 1 || !W.allFinite()) {
    return misskappa::unexpected(Error::dimension_mismatch);
  }
  if ((W - W.transpose()).cwiseAbs().maxCoeff() > 1e-10) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  Eigen::SelfAdjointEigenSolver<RealMat> es(0.5 * (W + W.transpose()));
  if (es.info() != Eigen::Success) return misskappa::unexpected(Error::numerical_error);
  if (es.eigenvalues().minCoeff() < -1e-8) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  return {};
}

// Brennan-Prediger uniform chance disagreement for vector-valued ratings with
// feature metric W: D_e^BP = sum_a W_aa * (range_a^2 / 6), the W-weighted sum of
// per-feature continuous-uniform baselines over the observed score range. This
// is the vector generalisation of the scalar continuous BP baseline (under the
// independent-uniform chance model the off-diagonal W terms drop out).
// `ranges`, when non-empty, fixes the per-feature score width (hi_a - lo_a)
// instead of reading it from the observed data, mirroring the scalar `values`
// hook: it makes the chance denominator deterministic so the Brennan-Prediger
// SE is exact and the estimand no longer drifts with the sample range. It must
// be length `features` (one width per feature); an empty vector keeps the
// observed-range behaviour.
Result<double> vector_bp_chance(RealMatView ratings, int features, const RealMat& W,
                                const RealVec& ranges) {
  const bool user_ranges = ranges.size() == static_cast<Eigen::Index>(features);
  if (ranges.size() != 0 && !user_ranges) {
    return misskappa::unexpected(Error::dimension_mismatch);
  }
  const int cols = static_cast<int>(ratings.cols());
  const int R = cols / features;
  double chance = 0.0;
  for (int a = 0; a < features; ++a) {
    double width;
    if (user_ranges) {
      width = ranges(a);
      if (!std::isfinite(width) || width <= 0.0) {
        return misskappa::unexpected(Error::invalid_argument);
      }
    } else {
      double lo = std::numeric_limits<double>::infinity();
      double hi = -std::numeric_limits<double>::infinity();
      for (int r = 0; r < R; ++r) {
        const int col = vector_col(r, a, features);
        for (int i = 0; i < ratings.rows(); ++i) {
          const double v = ratings(i, col);
          if (!std::isfinite(v)) continue;
          lo = std::min(lo, v);
          hi = std::max(hi, v);
        }
      }
      if (!std::isfinite(lo) || !std::isfinite(hi)) {
        return misskappa::unexpected(Error::invalid_argument);
      }
      width = hi - lo;
    }
    chance += W(a, a) * width * width / 6.0;
  }
  if (!std::isfinite(chance) || chance <= 1e-12) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  return chance;
}

GradientResult vector_quadratic_gradient(
    const RealVec& mu, const RealMat& sigma, int R, int features, const RealMat& W,
    double bp_chance) {
  const int q0 = static_cast<int>(mu.size());
  const int pstar = vech_size(q0);
  const RealMat B = vector_design(R, features, W, 0);
  const RealMat T = vector_design(R, features, W, 1);
  const RealMat Gm = vector_design(R, features, W, 2);

  const double tB = (B.array() * sigma.array()).sum();
  const double tT = (T.array() * sigma.array()).sum();
  const double tG = mu.dot(Gm * mu);
  const double dC = static_cast<double>(R - 1) * tT + static_cast<double>(R) * tG;
  const double nC = tB - tT;
  const double dF = static_cast<double>(R - 1) * (tT + tG);
  const double nF = tB - tT - tG;
  const double d_obs =
      (2.0 / static_cast<double>(R - 1)) * (tT + tG - tB / static_cast<double>(R));

  RealMat At = RealMat::Zero(3, q0 + pstar);
  At.block(0, q0, 1, pstar) = vech_gradient_from_matrix(B).transpose();
  At.block(1, q0, 1, pstar) = vech_gradient_from_matrix(T).transpose();
  At.block(2, 0, 1, q0) = (2.0 * (Gm * mu)).transpose();

  RealMat Jh = RealMat::Zero(3, 3);
  Jh(0, 0) = 1.0 / dC;
  Jh(0, 1) = -1.0 / dC - nC * static_cast<double>(R - 1) / (dC * dC);
  Jh(0, 2) = -nC * static_cast<double>(R) / (dC * dC);
  const double f_t23 = -1.0 / dF - nF * static_cast<double>(R - 1) / (dF * dF);
  Jh(1, 0) = 1.0 / dF;
  Jh(1, 1) = f_t23;
  Jh(1, 2) = f_t23;
  const double bp_factor = -2.0 / (bp_chance * static_cast<double>(R - 1));
  Jh(2, 0) = -bp_factor / static_cast<double>(R);
  Jh(2, 1) = bp_factor;
  Jh(2, 2) = bp_factor;

  GradientResult out;
  out.estimates = RealVec(3);
  out.estimates(0) = nC / dC;
  out.estimates(1) = nF / dF;
  out.estimates(2) = 1.0 - d_obs / bp_chance;
  out.gradient = Jh * At;
  return out;
}

}  // namespace

Result<NormalFimlEstimation> estimate_alpha_normal_fiml(
    RealMatView ratings, EmOptions opts) {
  return estimate_normal_fiml_impl(ratings, opts, [](const RealVec& mu, const RealMat& sigma) {
    return alpha_gradient(mu, sigma);
  });
}

Result<NormalFimlEstimation> estimate_quadratic_normal_fiml(
    RealMatView ratings, const RealVec& values, EmOptions opts) {
  auto bp_chance = quadratic_bp_chance(ratings, values);
  if (!bp_chance) return misskappa::unexpected(bp_chance.error());
  return estimate_normal_fiml_impl(
      ratings, opts, [chance = *bp_chance](const RealVec& mu, const RealMat& sigma) {
        return scalar_quadratic_gradient(mu, sigma, chance);
      });
}

Result<NormalFimlEstimation> estimate_vector_quadratic_normal_fiml(
    RealMatView ratings, int features, RealMatView W, EmOptions opts,
    const RealVec& ranges) {
  if (opts.covariance != NormalCovarianceModel::saturated &&
      opts.covariance != NormalCovarianceModel::variance_components) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  const int cols = static_cast<int>(ratings.cols());
  if (features < 1 || cols < 1 || cols % features != 0) {
    return misskappa::unexpected(Error::dimension_mismatch);
  }
  const int R = cols / features;
  if (R < 2) return misskappa::unexpected(Error::invalid_argument);
  if (W.rows() != features || W.cols() != features) {
    return misskappa::unexpected(Error::dimension_mismatch);
  }
  auto valid = validate_feature_metric(W);
  if (!valid) return misskappa::unexpected(valid.error());
  RealMat Wsym = 0.5 * (W + W.transpose());
  auto bp_chance = vector_bp_chance(ratings, features, Wsym, ranges);
  if (!bp_chance) return misskappa::unexpected(bp_chance.error());
  // The block compound-symmetry covariance Sigma = J_R (x) Psi + I_R (x) E is
  // selected by the feature block size; the saturated route keeps block_size 1.
  if (opts.covariance == NormalCovarianceModel::variance_components) {
    opts.block_size = features;
  }
  return estimate_normal_fiml_impl(
      ratings, opts,
      [R, features, Wsym, chance = *bp_chance](const RealVec& mu, const RealMat& sigma) {
        return vector_quadratic_gradient(mu, sigma, R, features, Wsym, chance);
      });
}

}  // namespace misskappa
