// Rcpp glue between R and the standalone C++17 misskappa library.
//
// The library is compiled with -fno-exceptions and returns Result<T> by value
// for fallible APIs. This TU runs with exceptions enabled (Rcpp requires it)
// and translates Result<T> errors into Rcpp::stop. The library itself never
// throws, so the no-exceptions / exceptions boundary is safe.

#include <Rcpp.h>
#include <RcppEigen.h>

#include <cmath>
#include <limits>

#include "misskappa/estimate.hpp"
#include "misskappa/diagnostics.hpp"
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
  Eigen::VectorXd values;
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

  return {std::move(mapped), std::move(W), std::move(v), C};
}

Rcpp::List estimation_to_list(const misskappa::Estimation& e) {
  Rcpp::NumericVector est(e.estimates.size());
  for (Eigen::Index i = 0; i < e.estimates.size(); ++i) est[i] = e.estimates(i);
  Rcpp::NumericMatrix vcov(e.vcov.rows(), e.vcov.cols());
  for (Eigen::Index i = 0; i < e.vcov.rows(); ++i) {
    for (Eigen::Index j = 0; j < e.vcov.cols(); ++j) vcov(i, j) = e.vcov(i, j);
  }
  // psi is the per-subject influence-function matrix (n x K). Estimators
  // that do not expose IFs leave it empty (0 x 0); that empty shape
  // round-trips to R as a 0 x 0 numeric matrix, which the R-side
  // accessor reads as "no IF available".
  Rcpp::NumericMatrix psi(e.psi.rows(), e.psi.cols());
  for (Eigen::Index i = 0; i < e.psi.rows(); ++i) {
    for (Eigen::Index j = 0; j < e.psi.cols(); ++j) psi(i, j) = e.psi(i, j);
  }
  return Rcpp::List::create(
      Rcpp::Named("estimates") = est,
      Rcpp::Named("vcov") = vcov,
      Rcpp::Named("psi") = psi);
}

misskappa::EmOptions parse_em_options(Rcpp::List em_options) {
  misskappa::EmOptions opts;
  if (em_options.containsElementNamed("tol"))
    opts.tol = Rcpp::as<double>(em_options["tol"]);
  if (em_options.containsElementNamed("max_iter"))
    opts.max_iter = Rcpp::as<int>(em_options["max_iter"]);
  if (em_options.containsElementNamed("prune_tol"))
    opts.prune_tol = Rcpp::as<double>(em_options["prune_tol"]);
  if (em_options.containsElementNamed("start_alpha"))
    opts.start_alpha = Rcpp::as<double>(em_options["start_alpha"]);
  if (em_options.containsElementNamed("info_rcond"))
    opts.info_rcond = Rcpp::as<double>(em_options["info_rcond"]);
  return opts;
}

Rcpp::NumericVector to_numeric_vector(const Eigen::VectorXd& x) {
  Rcpp::NumericVector out(x.size());
  for (Eigen::Index i = 0; i < x.size(); ++i) out[i] = x(i);
  return out;
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
    misskappa::EmOptions opts = parse_em_options(em_options);
    r = misskappa::estimate_fiml(in.ratings_indexed, in.weights, opts);
  } else {
    Rcpp::stop("Unknown method: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_alpha_raw(
    const Rcpp::IntegerMatrix& x,
    std::string method,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options) {
  PreparedInputs in = prepare_inputs(x, "identity", values);

  misskappa::Result<misskappa::Estimation> r;
  if (method == "available") {
    r = misskappa::estimate_alpha_available(in.ratings_indexed, in.values);
  } else if (method == "fiml") {
    misskappa::EmOptions opts = parse_em_options(em_options);
    r = misskappa::estimate_alpha_fiml(in.ratings_indexed, in.values, opts);
  } else {
    Rcpp::stop("Unknown method: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_alpha_available_continuous(const Rcpp::NumericMatrix& x) {
  Eigen::MatrixXd ratings(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) ratings(i, j) = x(i, j);
  }

  auto r = misskappa::estimate_alpha_available_continuous(ratings);
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_fiml_counts(
    const Rcpp::IntegerMatrix& x,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    int r_total,
    Rcpp::List em_options) {
  const int C = x.ncol();
  if (C < 1) Rcpp::stop("Counts matrix must have at least one category column.");

  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mapped(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) mapped(i, j) = x(i, j);
  }

  Eigen::VectorXd v(C);
  if (values.isNotNull()) {
    Rcpp::NumericVector vv(values);
    if (vv.size() != C) Rcpp::stop("Length of 'values' must equal the number of category columns.");
    for (int i = 0; i < C; ++i) v(i) = vv[i];
  } else {
    for (int i = 0; i < C; ++i) v(i) = static_cast<double>(i + 1);
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

  misskappa::EmOptions opts = parse_em_options(em_options);

  auto r = misskappa::estimate_fiml_counts(mapped, W, r_total, opts);
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_fiml_louis_spectrum(
    const Rcpp::IntegerMatrix& x,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values,
    Rcpp::List em_options) {
  PreparedInputs in = prepare_inputs(x, weight_type, values);
  misskappa::EmOptions opts = parse_em_options(em_options);
  auto r = unwrap(misskappa::diagnose_fiml_louis(in.ratings_indexed, in.weights, opts));
  return Rcpp::List::create(
      Rcpp::Named("eigenvalues") = to_numeric_vector(r.eigenvalues),
      Rcpp::Named("gradient_projection") = to_numeric_vector(r.gradient_projection),
      Rcpp::Named("variance_contribution") = to_numeric_vector(r.variance_contribution),
      Rcpp::Named("variance") = r.variance,
      Rcpp::Named("lambda_max") = r.lambda_max,
      Rcpp::Named("threshold") = r.threshold,
      Rcpp::Named("retained_rank") = r.retained_rank,
      Rcpp::Named("kappa_conger") = r.kappa_conger,
      Rcpp::Named("n_subjects") = static_cast<double>(r.n_subjects),
      Rcpp::Named("n_patterns") = static_cast<double>(r.n_patterns),
      Rcpp::Named("C") = r.c,
      Rcpp::Named("R") = r.R);
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_quadratic(
    const Rcpp::NumericMatrix& x,
    Rcpp::NumericVector values) {
  Eigen::MatrixXd ratings(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) ratings(i, j) = x(i, j);
  }
  Eigen::VectorXd v(values.size());
  for (int i = 0; i < values.size(); ++i) v(i) = values[i];

  auto r = misskappa::estimate_quadratic(ratings, v);
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_counts(
    const Rcpp::IntegerMatrix& x,
    std::string weight_type,
    Rcpp::Nullable<Rcpp::NumericVector> values) {
  const int C = x.ncol();
  if (C < 1) Rcpp::stop("Counts matrix must have at least one category column.");

  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mapped(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) mapped(i, j) = x(i, j);
  }

  Eigen::VectorXd v(C);
  if (values.isNotNull()) {
    Rcpp::NumericVector vv(values);
    if (vv.size() != C) Rcpp::stop("Length of 'values' must equal the number of category columns.");
    for (int i = 0; i < C; ++i) v(i) = vv[i];
  } else {
    for (int i = 0; i < C; ++i) v(i) = static_cast<double>(i + 1);
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

  auto r = misskappa::estimate_available_counts(mapped, W);
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_continuous(
    const Rcpp::NumericMatrix& x,
    std::string method,
    std::string weight_type) {
  Eigen::MatrixXd ratings(x.nrow(), x.ncol());
  double mn = std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) {
      const double v = x(i, j);
      ratings(i, j) = v;
      if (std::isfinite(v)) {
        if (v < mn) mn = v;
        if (v > mx) mx = v;
      }
    }
  }
  if (!std::isfinite(mn)) Rcpp::stop("All ratings are missing.");

  misskappa::loss::ContinuousLoss loss;
  if (weight_type == "identity" || weight_type == "unweighted") {
    loss = unwrap(misskappa::loss::identity_loss());
  } else if (weight_type == "linear") {
    loss = unwrap(misskappa::loss::linear_loss(mn, mx));
  } else if (weight_type == "quadratic") {
    loss = unwrap(misskappa::loss::quadratic_loss(mn, mx));
  } else if (weight_type == "radical") {
    loss = unwrap(misskappa::loss::radical_loss(mn, mx));
  } else if (weight_type == "ratio") {
    loss = unwrap(misskappa::loss::ratio_loss(mn, mx));
  } else {
    Rcpp::stop("Unknown weight type for continuous data: " + weight_type);
  }

  misskappa::Result<misskappa::Estimation> r;
  if (method == "available") {
    r = misskappa::estimate_available_continuous(ratings, loss);
  } else if (method == "ipw") {
    r = misskappa::estimate_ipw_continuous(ratings, loss);
  } else if (method == "gwet") {
    r = misskappa::estimate_gwet_continuous(ratings, loss);
  } else {
    Rcpp::stop("Unknown method for continuous data: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_vector(
    const Rcpp::NumericMatrix& x,
    int features,
    std::string method,
    std::string loss_type,
    Rcpp::NumericVector feature_weights) {
  Eigen::MatrixXd ratings(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) ratings(i, j) = x(i, j);
  }

  Eigen::VectorXd weights(feature_weights.size());
  for (int i = 0; i < feature_weights.size(); ++i) weights(i) = feature_weights[i];

  misskappa::loss::ComponentSeparableLoss loss;
  if (loss_type == "hamming") {
    loss = unwrap(misskappa::loss::hamming_vector_loss(weights));
  } else if (loss_type == "absolute") {
    loss = unwrap(misskappa::loss::absolute_vector_loss(weights));
  } else if (loss_type == "squared") {
    loss = unwrap(misskappa::loss::squared_vector_loss(weights));
  } else if (loss_type == "rms") {
    loss = unwrap(misskappa::loss::rms_vector_loss(weights));
  } else {
    Rcpp::stop("Unknown vector loss type: " + loss_type);
  }

  misskappa::Result<misskappa::Estimation> r;
  if (method == "pairwise") {
    r = misskappa::estimate_pairwise_vector(ratings, features, loss);
  } else if (method == "ipw") {
    r = misskappa::estimate_ipw_vector(ratings, features, loss);
  } else {
    Rcpp::stop("Unknown vector method: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_gwise_categorical(
    const Rcpp::IntegerMatrix& x,
    std::string distance_type,
    std::string method,
    int g,
    int max_chance_tuples,
    Rcpp::List em_options) {
  // The complete-data estimator requires every entry observed; the IPW and
  // FIML estimators accept missing entries (NA / na_code).
  const bool allow_missing = (method == "ipw" || method == "fiml");
  int C = 0;
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> ratings(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) {
      const int v = x(i, j);
      if (v == NA_INTEGER || v == kPkgNaCode) {
        if (!allow_missing) {
          Rcpp::stop("Complete-data g-wise estimator requires complete ratings; "
                     "use estimator = \"ipw\" or \"cat_fiml\" for missing data.");
        }
        ratings(i, j) = kPkgNaCode;
        continue;
      }
      if (v < 0) Rcpp::stop("G-wise categorical ratings must be non-negative codes.");
      ratings(i, j) = v;
      if (v + 1 > C) C = v + 1;
    }
  }
  if (C < 1) Rcpp::stop("All ratings are missing.");

  misskappa::loss::GwiseCategoricalDistance distance;
  if (distance_type == "nominal") {
    distance = unwrap(misskappa::loss::frechet_nominal_distance(C));
  } else if (distance_type == "hubert") {
    distance = unwrap(misskappa::loss::hubert_categorical_distance(C));
  } else {
    Rcpp::stop("Unknown categorical g-wise distance: " + distance_type);
  }

  misskappa::GwiseOptions opts;
  opts.g = g;
  opts.max_chance_tuples = max_chance_tuples;

  misskappa::Result<misskappa::Estimation> r;
  if (method == "complete") {
    r = misskappa::estimate_gwise(ratings, distance, opts);
  } else if (method == "ipw") {
    r = misskappa::estimate_ipw_gwise(ratings, distance, opts);
  } else if (method == "fiml") {
    misskappa::EmOptions em = parse_em_options(em_options);
    r = misskappa::estimate_fiml_gwise(ratings, distance, em, opts);
  } else {
    Rcpp::stop("Unknown categorical g-wise method: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}

// [[Rcpp::export]]
Rcpp::List rcpp_kappa_gwise_continuous(
    const Rcpp::NumericMatrix& x,
    std::string distance_type,
    std::string method,
    int g,
    int max_chance_tuples) {
  Eigen::MatrixXd ratings(x.nrow(), x.ncol());
  for (int i = 0; i < x.nrow(); ++i) {
    for (int j = 0; j < x.ncol(); ++j) ratings(i, j) = x(i, j);
  }

  misskappa::loss::GwiseContinuousDistance distance;
  if (distance_type == "absolute") {
    distance = unwrap(misskappa::loss::frechet_absolute_distance());
  } else if (distance_type == "quadratic") {
    distance = unwrap(misskappa::loss::frechet_quadratic_distance());
  } else if (distance_type == "hubert") {
    distance = unwrap(misskappa::loss::hubert_continuous_distance());
  } else {
    Rcpp::stop("Unknown continuous g-wise distance: " + distance_type);
  }

  misskappa::GwiseOptions opts;
  opts.g = g;
  opts.max_chance_tuples = max_chance_tuples;

  misskappa::Result<misskappa::Estimation> r;
  if (method == "complete") {
    r = misskappa::estimate_gwise_continuous(ratings, distance, opts);
  } else if (method == "ipw") {
    r = misskappa::estimate_ipw_gwise_continuous(ratings, distance, opts);
  } else {
    Rcpp::stop("Unknown continuous g-wise method: " + method);
  }
  return estimation_to_list(unwrap(std::move(r)));
}
