#include "doctest.h"

#include "detail_psd_inverse.hpp"
#include "misskappa/types.hpp"

#include <Eigen/QR>
#include <cmath>

using misskappa::RealMat;
using misskappa::RealVec;
namespace detail = misskappa::detail;

namespace {

constexpr double tol = 1e-9;

// A fixed orthonormal basis Q (columns q0, q1, q2) so the eigenstructure of the
// constructed information matrices is known exactly and platform-independent.
RealMat fixed_orthonormal() {
  RealMat a(3, 3);
  a << 1.0, 0.3, -0.7,
       0.2, 1.0, 0.4,
       -0.5, 0.6, 1.0;
  Eigen::HouseholderQR<RealMat> qr(a);
  return RealMat(qr.householderQ());
}

RealMat spectral(const RealMat& q, double l0, double l1, double l2) {
  RealVec lam(3);
  lam << l0, l1, l2;
  return q * lam.asDiagonal() * q.transpose();
}

}  // namespace

TEST_CASE("pseudo_inverse_psd: full-rank matches the ordinary inverse") {
  const RealMat q = fixed_orthonormal();
  const RealMat m = spectral(q, 4.0, 2.0, 1.0);  // SPD, well conditioned

  const RealMat pinv = detail::pseudo_inverse_psd(m, 1e-12);
  CHECK((pinv - m.inverse()).cwiseAbs().maxCoeff() < tol);
  // Analytic inverse: Q diag(1/4, 1/2, 1) Q^T.
  const RealMat expected = spectral(q, 0.25, 0.5, 1.0);
  CHECK((pinv - expected).cwiseAbs().maxCoeff() < tol);
  CHECK((pinv - pinv.transpose()).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("pseudo_inverse_psd: rank-deficient gives the Moore-Penrose inverse") {
  const RealMat q = fixed_orthonormal();
  const RealMat m = spectral(q, 4.0, 2.0, 0.0);  // rank 2, null space = span(q2)

  const RealMat pinv = detail::pseudo_inverse_psd(m, 1e-8);
  // The truncated direction is dropped, not inverted: Q diag(1/4, 1/2, 0) Q^T.
  const RealMat expected = spectral(q, 0.25, 0.5, 0.0);
  CHECK((pinv - expected).cwiseAbs().maxCoeff() < tol);

  // Penrose conditions for the pseudo-inverse.
  CHECK((m * pinv * m - m).cwiseAbs().maxCoeff() < tol);
  CHECK((pinv * m * pinv - pinv).cwiseAbs().maxCoeff() < tol);
  // M * M^+ is the orthogonal projector onto range(M) = I - q2 q2^T.
  const RealVec q2 = q.col(2);
  const RealMat proj = RealMat::Identity(3, 3) - q2 * q2.transpose();
  CHECK((m * pinv - proj).cwiseAbs().maxCoeff() < tol);
}

TEST_CASE("reduced_gradient_null_fraction: 0 in range, 1 in null, mix in between") {
  const RealMat q = fixed_orthonormal();
  const RealMat info = spectral(q, 4.0, 2.0, 0.0);  // null space = span(q2)

  // Each row is a unit "gradient" direction whose null-space content is known.
  RealMat jac(3, 3);
  jac.row(0) = q.col(0).transpose();                              // in range -> 0
  jac.row(1) = q.col(2).transpose();                              // in null  -> 1
  jac.row(2) = ((q.col(0) + q.col(2)).normalized()).transpose();  // 50/50 -> 1/sqrt2

  auto frac = detail::reduced_gradient_null_fraction(info, jac, 1e-8);
  REQUIRE(frac.has_value());
  CHECK(std::abs((*frac)(0) - 0.0) < tol);
  CHECK(std::abs((*frac)(1) - 1.0) < tol);
  CHECK(std::abs((*frac)(2) - std::sqrt(0.5)) < tol);
}

TEST_CASE("reduced_gradient_null_fraction: full-rank information reports no leak") {
  const RealMat q = fixed_orthonormal();
  const RealMat info = spectral(q, 4.0, 2.0, 1.0);  // full rank

  RealMat jac(2, 3);
  jac.row(0) = q.col(0).transpose();
  jac.row(1) = q.col(2).transpose();

  auto frac = detail::reduced_gradient_null_fraction(info, jac, 1e-8);
  REQUIRE(frac.has_value());
  CHECK((*frac).maxCoeff() < tol);  // nothing truncated -> zero everywhere
}
