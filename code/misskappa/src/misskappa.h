#ifndef MISSKAPPA_HPP
#define MISSKAPPA_HPP

#include "emdiscrete.h" // Our kappa library DEPENDS on the EM library

namespace misskappa {
// --- Core Kappa-specific Types ---
using uvec = emdiscrete::uvec;

template <typename T>
using Result = emdiscrete::Result<T>;

// Final estimation result struct returned by all kappa methods
struct Estimation {
  arma::vec estimates;
  arma::mat vcov;
};

// --- Loss/Weight Matrix Generation (for Categorical Data) ---
namespace loss {
Result<arma::mat> identity_weights(int c);
Result<arma::mat> linear_weights(int c, const arma::vec& v);
Result<arma::mat> quadratic_weights(int c, const arma::vec& v);
Result<arma::mat> ordinal_weights(int c);
Result<arma::mat> radical_weights(int c, const arma::vec& v);
Result<arma::mat> ratio_weights(int c, const arma::vec& v);
Result<arma::mat> circular_weights(int c, const arma::vec& v);
Result<arma::mat> bipolar_weights(int c, const arma::vec& v);

// --- NEW: Loss Function Generation (for Continuous Data) ---
using LossFunction = std::function<double(double, double)>;
Result<LossFunction> create_identity_loss();
Result<LossFunction> create_linear_loss(double min_val, double max_val);
Result<LossFunction> create_quadratic_loss(double min_val, double max_val);
Result<LossFunction> create_radical_loss(double min_val, double max_val);
Result<LossFunction> create_ratio_loss(double min_val, double max_val);
}

// --- Public API for Kappa Estimators ---
namespace kappaml {
Result<Estimation> kappa(const emdiscrete::EM_Result&, const arma::mat& loss_matrix);
Result<Estimation> kappa_counts(const emdiscrete::EM_Result&, const arma::mat& loss_matrix);
}
namespace kappanp {
Result<Estimation> kappa(const arma::imat&, const arma::mat& loss_matrix, bool use_ipw, bool use_gwet);
Result<Estimation> kappa_continuous(const arma::mat&, const loss::LossFunction& loss_func, bool use_ipw, bool use_gwet);
Result<Estimation> kappa_counts(const arma::umat& counts, const arma::mat& loss_matrix);

}
namespace kappaqp {
Result<Estimation> kappa(const arma::mat& x, const arma::vec& values);
Result<Estimation> kappa_counts(const arma::mat& x, const arma::vec& values, int R);
}
}
#endif // MISSKAPPA_HPP
