// Phase-timing accumulators for the instrumented estimate_fiml copy.
#ifndef CATFIML_BENCH_PROF_HPP
#define CATFIML_BENCH_PROF_HPP

#include <chrono>

namespace bench_prof {

struct Accum {
  double pre = 0, em = 0, louis = 0, pinv = 0, theta_vcov = 0,
         kappa_map = 0, null_frac = 0, rest = 0;
  long em_iters = 0;
  long n_active = 0;   // cells in the EM state vector (currently C^R)
  long n_final = 0;    // pruned support size
  long n_groups = 0;
  void reset() { *this = Accum{}; }
};

inline Accum& acc() {
  static Accum a;
  return a;
}

struct Timer {
  double& slot;
  std::chrono::steady_clock::time_point t0;
  explicit Timer(double& s) : slot(s), t0(std::chrono::steady_clock::now()) {}
  ~Timer() {
    slot += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
  }
};

}  // namespace bench_prof

#endif
