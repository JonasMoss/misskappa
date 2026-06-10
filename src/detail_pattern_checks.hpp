#ifndef MISSKAPPA_SRC_DETAIL_PATTERN_CHECKS_HPP
#define MISSKAPPA_SRC_DETAIL_PATTERN_CHECKS_HPP

#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace misskappa::detail {

inline Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
observed_mask(IntMatView ratings) {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mask(ratings.rows(), ratings.cols());
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      mask(i, j) = (ratings(i, j) == na_code) ? 0 : 1;
    }
  }
  return mask;
}

inline Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>
finite_mask(RealMatView ratings) {
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mask(ratings.rows(), ratings.cols());
  for (Eigen::Index i = 0; i < ratings.rows(); ++i) {
    for (Eigen::Index j = 0; j < ratings.cols(); ++j) {
      mask(i, j) = std::isfinite(ratings(i, j)) ? 1 : 0;
    }
  }
  return mask;
}

inline Result<void> require_complete_pair_observation(
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mask) {
  const Eigen::Index n = mask.rows();
  const Eigen::Index R = mask.cols();
  if (n < 1 || R < 2) return misskappa::unexpected(Error::invalid_argument);

  for (Eigen::Index j = 0; j < R; ++j) {
    bool seen = false;
    for (Eigen::Index i = 0; i < n; ++i) {
      if (mask(i, j) != 0) {
        seen = true;
        break;
      }
    }
    if (!seen) return misskappa::unexpected(Error::not_identified);
  }

  for (Eigen::Index j = 0; j < R - 1; ++j) {
    for (Eigen::Index k = j + 1; k < R; ++k) {
      bool jointly_seen = false;
      for (Eigen::Index i = 0; i < n; ++i) {
        if (mask(i, j) != 0 && mask(i, k) != 0) {
          jointly_seen = true;
          break;
        }
      }
      if (!jointly_seen) return misskappa::unexpected(Error::not_identified);
    }
  }
  return {};
}

inline Result<void> require_complete_tuple_observation(
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& mask, int tuple_size) {
  const Eigen::Index n = mask.rows();
  const Eigen::Index R = mask.cols();
  if (n < 1 || R < 2 || tuple_size < 2 || tuple_size > R) {
    return misskappa::unexpected(Error::invalid_argument);
  }
  if (tuple_size == 2) return require_complete_pair_observation(mask);

  std::vector<Eigen::Index> tuple(static_cast<std::size_t>(tuple_size));
  for (int g = 0; g < tuple_size; ++g) tuple[static_cast<std::size_t>(g)] = g;

  while (true) {
    bool jointly_seen = false;
    for (Eigen::Index i = 0; i < n; ++i) {
      bool row_complete = true;
      for (Eigen::Index j : tuple) {
        if (mask(i, j) == 0) {
          row_complete = false;
          break;
        }
      }
      if (row_complete) {
        jointly_seen = true;
        break;
      }
    }
    if (!jointly_seen) return misskappa::unexpected(Error::not_identified);

    int pos = tuple_size - 1;
    while (pos >= 0) {
      const Eigen::Index limit = R - tuple_size + pos;
      if (tuple[static_cast<std::size_t>(pos)] < limit) break;
      --pos;
    }
    if (pos < 0) break;

    ++tuple[static_cast<std::size_t>(pos)];
    for (int g = pos + 1; g < tuple_size; ++g) {
      tuple[static_cast<std::size_t>(g)] = tuple[static_cast<std::size_t>(g - 1)] + 1;
    }
  }
  return {};
}

}  // namespace misskappa::detail

#endif  // MISSKAPPA_SRC_DETAIL_PATTERN_CHECKS_HPP
