// Rcpp glue between R and the standalone C++23 misskappa library.
//
// The library is compiled with -fno-exceptions and returns Result<T> by value
// for fallible APIs. This TU runs with exceptions enabled (Rcpp requires it)
// and translates Result<T> errors into Rcpp::stop. The library itself never
// throws, so the no-exceptions / exceptions boundary is safe.

#include <Rcpp.h>
#include <RcppEigen.h>

#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"

// [[Rcpp::depends(RcppEigen)]]

namespace {

constexpr int kPkgNaCode = -1;

const char* error_to_message(misskappa::Error e) {
  using E = misskappa::Error;
  switch (e) {
    case E::invalid_argument:   return "Invalid argument.";
    case E::dimension_mismatch: return "Dimension mismatch.";
    case E::singular_weight:    return "Rater has zero observations; IPW/Gwet weight is singular.";
    case E::numerical_error:    return "Numerical error.";
    case E::not_supported:      return "Operation not supported.";
    case E::not_converged:      return "EM did not converge within max_iter iterations.";
  }
  return "Unknown error.";
}

template <typename T>
T unwrap(misskappa::Result<T>&& r) {
  if (!r.has_value()) Rcpp::stop(error_to_message(r.error()));
  return std::move(*r);
}

// Translate an R integer matrix to Eigen, mapping NA_INTEGER to na_code.
// Also remaps non-negative category values to a sorted-unique index range
// [0, C-1], returns the C x C agreement weight matrix built from the
// requested weighting and the original category values.
struct PreparedInputs {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> ratings_indexed;
  Eigen::MatrixXd weights;
  int C;
};

PreparedInputs prepare_inputs(
    const Rcpp::IntegerMatrix& x,
    const std::string& weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values) {
  // Collect unique non-NA category values.
  std::set<int> cat_set;
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) {
      const int v = x(i, j);
      if (v == NA_INTEGER) continue;
      cat_set.insert(v);
    }
  }
  if (cat_set.empty()) Rcpp::stop("All ratings are missing.");
  std::vector<int> cats(cat_set.begin(), cat_set.end());
  const int C = static_cast<int>(cats.size());

  std::map<int, int> cat_to_idx;
  for (int i = 0; i < C; ++i) cat_to_idx[cats[i]] = i;

  // Remap to [0, C-1] / na_code.
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mapped(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) {
      const int v = x(i, j);
      mapped(i, j) = (v == NA_INTEGER) ? kPkgNaCode : cat_to_idx.at(v);
    }
  }

  // Category value vector for metric weighting schemes. Default to the
  // category integers themselves (1, 2, ..., C if user passed 1-indexed,
  // 0..C-1 if user passed 0-indexed). irrCAC convention is 1-indexed.
  Eigen::VectorXd v(C);
  if (values.isNotNull()) {
    Rcpp::NumericVector vv(values);
    if (vv.size() != C) Rcpp::stop("Length of 'values' must equal the number of unique categories.");
    for (int i = 0; i < C; ++i) v(i) = vv[i];
  } else {
    for (int i = 0; i < C; ++i) v(i) = static_cast<double>(cats[i]);
  }

  Eigen::MatrixXd W;
  if (weight_type == "identity" || weight_type == "unweighted") {
    W = unwrap(misskappa::loss::identity_weights(C));
  } else if (weight_type == "linear") {
    W = unwrap(misskappa::loss::linear_weights(C, v));
  } else if (weight_type == "quadratic") {
    W = unwrap(misskappa::loss::quadratic_weights(C, v));
  } else if (weight_type == "ordinal") {
    W = unwrap(misskappa::loss::ordinal_weights(C));
  } else if (weight_type == "radical") {
    W = unwrap(misskappa::loss::radical_weights(C, v));
  } else if (weight_type == "ratio") {
    W = unwrap(misskappa::loss::ratio_weights(C, v));
  } else if (weight_type == "circular") {
    W = unwrap(misskappa::loss::circular_weights(C, v));
  } else if (weight_type == "bipolar") {
    W = unwrap(misskappa::loss::bipolar_weights(C, v));
  } else {
    Rcpp::stop("Unknown weight type: " + weight_type);
  }

  return {std::move(mapped), std::move(W), C};
}

Rcpp::List estimation_to_list(const misskappa::Estimation& e) {
  Rcpp::NumericVector est(e.estimates.size());
  for (Eigen::Index i = 0; i < e.estimates.size(); ++i) est[i] = e.estimates(i);
  Rcpp::NumericMatrix vcov(e.vcov.rows(), e.vcov.cols());
  for (Eigen::Index i = 0; i < e.vcov.rows(); ++i) {
    for (Eigen::Index j = 0; j < e.vcov.cols(); ++j) vcov(i, j) = e.vcov(i, j);
  }
  return Rcpp::List::create(
      Rcpp::Named("estimates") = est,
      Rcpp::Named("vcov") = vcov);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_raw(
    const Rcpp::IntegerMatrix& x,
    std::string method,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options) {
  PreparedInputs in = prepare_inputs(x, weight_type, values);

  misskappa::Result<misskappa::Estimation> r;
  if (method == "available") {
    r = misskappa::estimate_available(in.ratings_indexed, in.weights);
  } else if (method == "ipw") {
    r = misskappa::estimate_ipw(in.ratings_indexed, in.weights);
  } else if (method == "gwet") {
    r = misskappa::estimate_gwet(in.ratings_indexed, in.weights);
  } else if (method == "fiml") {
    misskappa::EmOptions opts;
    if (em_options.containsElementNamed("tol"))
      opts.tol = Rcpp::as<double>(em_options["tol"]);
    if (em_options.containsElementNamed("max_iter"))
      opts.max_iter = Rcpp::as<int>(em_options["max_iter"]);
    if (em_options.containsElementNamed("prune_tol"))
      opts.prune_tol = Rcpp::as<double>(em_options["prune_tol"]);
    if (em_options.containsElementNamed("start_alpha"))
      opts.start_alpha = Rcpp::as<double>(em_options["start_alpha"]);
    r = misskappa::estimate_fiml(in.ratings_indexed, in.weights, opts);
  } else {
    Rcpp::stop("Unknown method: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}
