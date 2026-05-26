#ifndef EMDISCRETE_HPP
#define EMDISCRETE_HPP

#include <RcppArmadillo.h>
#include <functional>
#include <string>
#include <vector>

#include "core/result.h"
#include "core/types.h"

namespace emdiscrete {

using uvec = misskappa_core::uvec;
constexpr int kNaInteger = misskappa_core::kNaInteger;
using Status = misskappa_core::Status;
template <typename T>
using Result = misskappa_core::Result<T>;

struct EM_Input {
  size_t n_total_patterns;
  int c;
  int r;
  std::string type;
  uvec pattern_indices_map;
  uvec group_n_subjects;
  uvec group_offsets;
  uvec group_n_completions;
  uvec completion_indices;
  arma::vec multivariate_hypergeom_coeffs;
};

struct EM_Options {
  double tol = 1e-8;
  int max_iter = 10000;
  double prune_tol = 1e-9;
  double start_alpha = 0.1;
};

struct EM_Result {
  arma::vec theta_hat;
  arma::mat var; // Covariance of theta_hat parameters
  uvec pattern_indices;
  std::string type;
  int c;
  int r;
  int iterations;
  bool converged;
  arma::uword n_subjects;
};

Result<EM_Input> preprocess_raw(const arma::imat& x, int c);
Result<EM_Input> preprocess_counts(const arma::umat& x, int r);
Result<EM_Result> run_em(const EM_Input& em_input, const EM_Options& options);

} // namespace emdiscrete
#endif // EMDISCRETE_HPP
