#ifndef MISSKAPPA_DIAGNOSTICS_HPP
#define MISSKAPPA_DIAGNOSTICS_HPP

#include "misskappa/estimate.hpp"
#include "misskappa/result.hpp"
#include "misskappa/types.hpp"

#include <cstddef>

namespace misskappa {

// Internal diagnostic surface used by repo experiments. Not part of the
// stable package API.
struct FimlLouisDiagnostic {
  RealVec eigenvalues;               // descending eigenvalues of I_obs.
  RealVec gradient_projection;       // q_k^T grad for Conger's kappa.
  RealVec variance_contribution;     // (q_k^T grad)^2 / lambda_k, zero if dropped.
  double variance = 0.0;
  double lambda_max = 0.0;
  double threshold = 0.0;
  double kappa_conger = 0.0;
  int retained_rank = 0;
  int c = 0;
  int R = 0;
  std::size_t n_subjects = 0;
  std::size_t n_patterns = 0;
};

Result<FimlLouisDiagnostic> diagnose_fiml_louis(
    IntMatView ratings, RealMatView weights, EmOptions opts);

}  // namespace misskappa

#endif  // MISSKAPPA_DIAGNOSTICS_HPP
