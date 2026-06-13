// Probe 39 — ten Hove, Jorgensen & van der Ark (2025) MLE-RE estimator core.
//
// Two-way crossed random-effects model for interrater reliability:
//   y_{sr} = mu + a_s + b_r + e_{sr},
//   a_s ~ N(0, s2s),  b_r ~ N(0, s2r),  e_{sr} ~ N(0, s2sr),
// one observation per OBSERVED (s, r) cell; missing cells simply absent
// (ignorable planned-missing design). This is lme4's
//   lmer(y ~ 1 + (1|subject) + (1|rater), REML = FALSE).
//
// Point estimate only (no variance) for now.
//
// The marginal covariance of the stacked observed vector is
//   V = s2sr * I_n + Z G Z',   Z = [Z_s | Z_r]  (n x (S+R)),
//   G = blkdiag(s2s I_S, s2r I_R).
// Woodbury / matrix-determinant-lemma collapse every evaluation to one
// Cholesky of the (S+R) x (S+R) matrix  M = G^{-1} + (1/s2sr) Z'Z:
//   a' V^{-1} b = (1/s2sr) a'b - (1/s2sr^2) (Z'a)' M^{-1} (Z'b),
//   log|V|      = n log s2sr + log|G| + log|M|.
// mu is profiled out by GLS.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Marginal log-likelihood of the model at given (positive) variance
// components, with mu profiled out. Returns -Inf if M is not PD.
// [[Rcpp::export]]
List icc_loglik(double s2s, double s2r, double s2sr,
                IntegerVector subj, IntegerVector rat, NumericVector y,
                int S, int R) {
  const int n = y.size();
  const int q = S + R;

  VectorXd ns = VectorXd::Zero(S), nr = VectorXd::Zero(R);  // cell counts
  VectorXd sy = VectorXd::Zero(S), ry = VectorXd::Zero(R);  // y sums
  MatrixXd A  = MatrixXd::Zero(S, R);                        // incidence
  double sumy = 0.0, sumy2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const int s = subj[i] - 1, r = rat[i] - 1;  // R is 1-based
    ns[s] += 1.0; nr[r] += 1.0;
    sy[s] += y[i]; ry[r] += y[i];
    A(s, r) += 1.0;
    sumy += y[i]; sumy2 += y[i] * y[i];
  }

  const double inv = 1.0 / s2sr;
  MatrixXd M = MatrixXd::Zero(q, q);
  for (int s = 0; s < S; ++s) M(s, s)         = 1.0 / s2s + ns[s] * inv;
  for (int r = 0; r < R; ++r) M(S + r, S + r) = 1.0 / s2r + nr[r] * inv;
  M.topRightCorner(S, R)   = A * inv;
  M.bottomLeftCorner(R, S) = A.transpose() * inv;

  VectorXd zy(q), zx(q);
  zy << sy, ry;   // Z' y
  zx << ns, nr;   // Z' 1  (X = column of ones)

  Eigen::LLT<MatrixXd> llt(M);
  if (llt.info() != Eigen::Success)
    return List::create(_["loglik"] = R_NegInf, _["mu"] = NA_REAL);

  const VectorXd Mzy = llt.solve(zy);
  const VectorXd Mzx = llt.solve(zx);
  const double inv2 = inv * inv;
  const double XtViX = inv * n     - inv2 * zx.dot(Mzx);
  const double XtViy = inv * sumy  - inv2 * zx.dot(Mzy);
  const double ytViy = inv * sumy2 - inv2 * zy.dot(Mzy);

  const double mu   = XtViy / XtViX;
  const double quad = ytViy - XtViy * XtViy / XtViX;  // (y-Xmu)'V^{-1}(y-Xmu)

  MatrixXd L = llt.matrixL();
  const double logdetM = 2.0 * L.diagonal().array().log().sum();
  const double logdetG = S * std::log(s2s) + R * std::log(s2r);
  const double logdetV = n * std::log(s2sr) + logdetG + logdetM;

  const double LOG2PI = std::log(2.0 * M_PI);
  const double ll = -0.5 * (n * LOG2PI + logdetV + quad);
  return List::create(_["loglik"] = ll, _["mu"] = mu);
}
