#include "misskappa.h"

namespace misskappa {
namespace loss {

// All weight functions return AGREEMENT matrices. The caller converts to loss if needed.
Result<arma::mat> identity_weights(int c) { return {Status::kOk, arma::eye<arma::mat>(c, c), ""}; }

Result<arma::mat> linear_weights(int c, const arma::vec& v) {
  if (v.n_elem != static_cast<arma::uword>(c)) return {Status::kError, std::nullopt, "Length of 'values' must equal 'c'."};
  arma::mat W(c, c);
  double v_min = v.min(), v_max = v.max();
  double range = std::abs(v_max - v_min);
  if (range < 1e-9) return identity_weights(c);
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) { W(i,j) = 1.0 - std::abs(v(i) - v(j)) / range; }
  return {Status::kOk, W, ""};
}

Result<arma::mat> quadratic_weights(int c, const arma::vec& v) {
  if (v.n_elem != static_cast<arma::uword>(c)) return {Status::kError, std::nullopt, "Length of 'values' must equal 'c'."};
  arma::mat W(c, c);
  double v_min = v.min(), v_max = v.max();
  double range_sq = std::pow(v_max - v_min, 2);
  if (range_sq < 1e-9) return identity_weights(c);
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) { W(i,j) = 1.0 - std::pow(v(i) - v(j), 2) / range_sq; }
  return {Status::kOk, W, ""};
}

Result<arma::mat> ordinal_weights(int c) {
  arma::mat W(c, c);
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) {
    double nkl = std::abs(i - j) + 1.0;
    W(i,j) = nkl * (nkl - 1.0) / 2.0;
  }
  double max_w = W.max();
  if (max_w < 1e-9) return identity_weights(c);
  return {Status::kOk, 1.0 - W / max_w, ""};
}

Result<arma::mat> radical_weights(int c, const arma::vec& v) {
  if (v.n_elem != static_cast<arma::uword>(c)) return {Status::kError, std::nullopt, "Length of 'values' must equal 'c'."};
  arma::mat W(c, c);
  double v_min = v.min(), v_max = v.max();
  double range_sqrt = std::sqrt(std::abs(v_max - v_min));
  if (range_sqrt < 1e-9) return identity_weights(c);
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) { W(i,j) = 1.0 - std::sqrt(std::abs(v(i) - v(j))) / range_sqrt; }
  return {Status::kOk, W, ""};
}

Result<arma::mat> ratio_weights(int c, const arma::vec& v) {
  if (v.n_elem != static_cast<arma::uword>(c)) return {Status::kError, std::nullopt, "Length of 'values' must equal 'c'."};
  arma::mat W(c, c);
  double v_min = v.min(), v_max = v.max();
  if (std::abs(v_max + v_min) < 1e-9) return identity_weights(c);
  double den_term = (v_max - v_min) / (v_max + v_min);
  double den_sq = std::pow(den_term, 2);
  if (den_sq < 1e-9) return identity_weights(c);
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) {
    if (std::abs(v(i) + v(j)) < 1e-9) W(i,j) = 1.0;
    else {
      double num_term = (v(i) - v(j)) / (v(i) + v(j));
      W(i,j) = 1.0 - std::pow(num_term, 2) / den_sq;
    }
  }
  return {Status::kOk, W, ""};
}

Result<arma::mat> circular_weights(int c, const arma::vec& v) {
  if (v.n_elem != static_cast<arma::uword>(c)) return {Status::kError, std::nullopt, "Length of 'values' must equal 'c'."};
  arma::mat W(c, c);
  double v_min = v.min(), v_max = v.max();
  double U = v_max - v_min + 1.0;
  if (U < 1e-9) return identity_weights(c);
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) { W(i,j) = std::pow(std::sin(M_PI * (v(i) - v(j)) / U), 2); }
  double max_w = W.max();
  if (max_w < 1e-9) return identity_weights(c);
  return {Status::kOk, 1.0 - W / max_w, ""};
}

Result<arma::mat> bipolar_weights(int c, const arma::vec& v) {
  if (v.n_elem != static_cast<arma::uword>(c)) return {Status::kError, std::nullopt, "Length of 'values' must equal 'c'."};
  arma::mat W(c, c, arma::fill::zeros);
  double v_min = v.min(), v_max = v.max();
  for(int i=0; i<c; ++i) for(int j=0; j<c; ++j) {
    if (i != j) {
      double den = ((v(i) + v(j)) - 2 * v_min) * (2 * v_max - (v(i) + v(j)));
      if (std::abs(den) > 1e-9) { W(i,j) = std::pow(v(i) - v(j), 2) / den; }
    }
  }
  double max_w = W.max();
  if (max_w < 1e-9) return identity_weights(c);
  return {Status::kOk, 1.0 - W / max_w, ""};
}

// --- NEW: Implementations for Continuous Loss Function Factories ---
Result<LossFunction> create_identity_loss() {
  return {Status::kOk, [](double v1, double v2) -> double {
    return (std::abs(v1 - v2) < 1e-9) ? 0.0 : 1.0;
  }, ""};
}

Result<LossFunction> create_linear_loss(double min_val, double max_val) {
  double range = std::abs(max_val - min_val);
  if (range < 1e-9) return create_identity_loss();
  return {Status::kOk, [range](double v1, double v2) -> double {
    return std::abs(v1 - v2) / range;
  }, ""};
}

Result<LossFunction> create_quadratic_loss(double min_val, double max_val) {
  double range_sq = std::pow(max_val - min_val, 2);
  if (range_sq < 1e-9) return create_identity_loss();
  return {Status::kOk, [range_sq](double v1, double v2) -> double {
    return std::pow(v1 - v2, 2) / range_sq;
  }, ""};
}

Result<LossFunction> create_radical_loss(double min_val, double max_val) {
  double range_sqrt = std::sqrt(std::abs(max_val - min_val));
  if (range_sqrt < 1e-9) return create_identity_loss();
  return {Status::kOk, [range_sqrt](double v1, double v2) -> double {
    return std::sqrt(std::abs(v1 - v2)) / range_sqrt;
  }, ""};
}

Result<LossFunction> create_ratio_loss(double min_val, double max_val) {
  if (std::abs(max_val + min_val) < 1e-9) return create_identity_loss();
  double den_term = (max_val - min_val) / (max_val + min_val);
  double den_sq = std::pow(den_term, 2);
  if (den_sq < 1e-9) return create_identity_loss();
  return {Status::kOk, [den_sq](double v1, double v2) -> double {
    if (std::abs(v1 + v2) < 1e-9) return 0.0;
    double num_term = (v1 - v2) / (v1 + v2);
    return std::pow(num_term, 2) / den_sq;
  }, ""};
}

} // namespace loss
} // namespace misskappa
