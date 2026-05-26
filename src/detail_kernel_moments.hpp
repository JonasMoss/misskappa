#ifndef MISSKAPPA_SRC_DETAIL_KERNEL_MOMENTS_HPP
#define MISSKAPPA_SRC_DETAIL_KERNEL_MOMENTS_HPP

#include "misskappa/types.hpp"

namespace misskappa::detail {

// V-statistic kernel accumulator: row / column / total sums over an n x n
// kernel table, plus the linear influence-function formula
// phi_i = (mean_j K(i, j) - psi) + (mean_j K(j, i) - psi). Shared between
// estimate_raw.cpp and estimate_continuous.cpp, which build kernels with
// identical structure.
struct KernelMoments {
  RealVec row_sum;
  RealVec col_sum;
  double total = 0.0;

  explicit KernelMoments(int n)
      : row_sum(RealVec::Zero(n)), col_sum(RealVec::Zero(n)) {}

  void add(int row, int col, double value) {
    row_sum(row) += value;
    col_sum(col) += value;
    total += value;
  }

  double mean(int n) const {
    return total / (static_cast<double>(n) * n);
  }

  RealVec influence(double psi, int n) const {
    const double inv_n = 1.0 / static_cast<double>(n);
    return ((row_sum.array() * inv_n - psi)
            + (col_sum.array() * inv_n - psi)).matrix();
  }
};

}  // namespace misskappa::detail

#endif  // MISSKAPPA_SRC_DETAIL_KERNEL_MOMENTS_HPP
