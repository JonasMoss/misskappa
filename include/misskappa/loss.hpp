#ifndef MISSKAPPA_LOSS_HPP
#define MISSKAPPA_LOSS_HPP

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

namespace misskappa::loss {

// Weight matrices for categorical agreement. Each returns a c x c matrix W
// with W(i, j) = 1 when category i agrees with category j and W(i, j) < 1
// for partial agreement, modulo the sign conventions documented per factory.
// All entries on the diagonal are 1; all matrices are symmetric.
//
// `c` is the number of categories. `v` is a category value vector (length c)
// for the metric-aware weightings (linear, quadratic, radical, ratio,
// circular, bipolar). The vector encodes the metric placement of categories
// on the line / on the circle.
Result<RealMat> identity_weights(int c);
Result<RealMat> linear_weights(int c, const RealVec& v);
Result<RealMat> quadratic_weights(int c, const RealVec& v);
Result<RealMat> ordinal_weights(int c);
Result<RealMat> radical_weights(int c, const RealVec& v);
Result<RealMat> ratio_weights(int c, const RealVec& v);
Result<RealMat> circular_weights(int c, const RealVec& v);
Result<RealMat> bipolar_weights(int c, const RealVec& v);

// Continuous loss kernels: agreement up to a metric on R. Returned as a small
// POD wrapper around a plain function pointer plus parameter pack, deliberately
// not std::function (allocates, exception-tainted, virtual dispatch).
struct ContinuousLoss {
  double min_val;
  double max_val;
  double (*compute)(double a, double b, double min_val, double max_val);
};

Result<ContinuousLoss> identity_loss();
Result<ContinuousLoss> linear_loss(double min_val, double max_val);
Result<ContinuousLoss> quadratic_loss(double min_val, double max_val);
Result<ContinuousLoss> radical_loss(double min_val, double max_val);
Result<ContinuousLoss> ratio_loss(double min_val, double max_val);

// Symmetric g-wise disagreement kernels for complete rectangular designs.
// These are deliberately small POD wrappers around plain function pointers.
// Categorical kernels expect category codes in [0, C-1]. Continuous kernels
// expect finite real values. The estimator supplies exactly `g` values.
struct GwiseCategoricalDistance {
  int C;
  double (*compute)(const int* values, int g, int C);
};

struct GwiseContinuousDistance {
  double (*compute)(const double* values, int g);
};

Result<GwiseCategoricalDistance> frechet_nominal_distance(int C);
Result<GwiseCategoricalDistance> hubert_categorical_distance(int C);

Result<GwiseContinuousDistance> frechet_absolute_distance();
Result<GwiseContinuousDistance> frechet_quadratic_distance();
Result<GwiseContinuousDistance> hubert_continuous_distance();

}  // namespace misskappa::loss

#endif  // MISSKAPPA_LOSS_HPP
