// Question: where does estimate_fiml (cat_fiml raw path) spend its time on
// study-29-shaped workloads, and is there headroom?
//
// Workload: C categories x R raters, n subjects, with either MCAR(30%) or a
// designed two-observed-raters pattern (heaviest completion fan-out: C^3
// completions per subject). Strict ML (flatten = 0), default EmOptions.
//
// Build/run: see run.sh in this folder. Compiled -O3 -g against the canonical
// src/ so perf can attribute internal phases.

#include "misskappa/estimate.hpp"
#include "prof.hpp"

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using Clock = std::chrono::steady_clock;

namespace {

// Simple correlated DGP: latent subject category + per-rater noise.
misskappa::IntMat simulate(int n, int R, int C, std::mt19937& rng) {
  misskappa::IntMat x(n, R);
  std::uniform_int_distribution<int> cat(0, C - 1);
  std::bernoulli_distribution agree(0.7);
  for (int i = 0; i < n; ++i) {
    const int truth = cat(rng);
    for (int j = 0; j < R; ++j) {
      x(i, j) = agree(rng) ? truth : cat(rng);
    }
  }
  return x;
}

void apply_mcar(misskappa::IntMat& x, double p, std::mt19937& rng) {
  std::bernoulli_distribution miss(p);
  for (int i = 0; i < x.rows(); ++i) {
    int dropped = 0;
    for (int j = 0; j < x.cols(); ++j) {
      if (dropped < x.cols() - 2 && miss(rng)) {
        x(i, j) = misskappa::na_code;
        ++dropped;
      }
    }
  }
}

// Each subject keeps exactly 2 raters, rotating so all pairs co-observe.
void apply_designed_random2(misskappa::IntMat& x, std::mt19937& rng) {
  const int R = static_cast<int>(x.cols());
  std::vector<std::pair<int, int>> pairs;
  for (int a = 0; a < R - 1; ++a)
    for (int b = a + 1; b < R; ++b) pairs.emplace_back(a, b);
  std::uniform_int_distribution<std::size_t> pick(0, pairs.size() - 1);
  for (int i = 0; i < x.rows(); ++i) {
    const auto [a, b] = pairs[i < static_cast<int>(pairs.size())
                                  ? static_cast<std::size_t>(i)  // guarantee coverage
                                  : pick(rng)];
    for (int j = 0; j < R; ++j) {
      if (j != a && j != b) x(i, j) = misskappa::na_code;
    }
  }
}

double run_case(const char* label, int n, int R, int C, bool designed,
                int reps, unsigned seed, double flatten = 0.0) {
  std::mt19937 rng(seed);
  const misskappa::RealMat W = misskappa::RealMat::Identity(C, C);
  // Match study-29's cat_fiml options.
  misskappa::EmOptions opts;
  opts.tol = 1e-7;
  opts.max_iter = 12000;
  opts.prune_tol = 1e-10;
  opts.flatten = flatten;

  bench_prof::acc().reset();
  double total_s = 0.0;
  double checksum = 0.0;
  double est_sum = 0.0;
  double se_sum = 0.0;
  int ok = 0;
  for (int r = 0; r < reps; ++r) {
    misskappa::IntMat x = simulate(n, R, C, rng);
    if (designed) {
      apply_designed_random2(x, rng);
    } else {
      apply_mcar(x, 0.3, rng);
    }
    const auto t0 = Clock::now();
    auto fit = misskappa::estimate_fiml(x, W, opts);
    const auto t1 = Clock::now();
    total_s += std::chrono::duration<double>(t1 - t0).count();
    if (fit) {
      ++ok;
      checksum += fit->estimates.sum() + fit->vcov.sum();
      est_sum += fit->estimates.sum();
      se_sum += fit->vcov.diagonal().cwiseSqrt().sum();
    }
  }
  const auto& a = bench_prof::acc();
  const double misc = a.rest - a.kappa_map - a.null_frac;
  std::printf("%-26s n=%4d R=%d C=%d reps=%d ok=%d mean=%8.1f ms est=%.10g se=%.10g\n",
              label, n, R, C, reps, ok, 1000.0 * total_s / reps, est_sum, se_sum);
  (void)checksum;
  std::printf("    state=%ld pruned=%ld groups=%ld em_iters/fit=%ld\n",
              a.n_active, a.n_final, a.n_groups, a.em_iters / reps);
  std::printf("    pre=%.2fs em=%.2fs louis=%.2fs pinv=%.2fs thvcov=%.2fs "
              "kmap=%.2fs nullfrac=%.2fs misc=%.2fs | sum=%.2fs total=%.2fs\n",
              a.pre, a.em, a.louis, a.pinv, a.theta_vcov, a.kappa_map,
              a.null_frac, misc,
              a.pre + a.em + a.louis + a.pinv + a.theta_vcov + a.kappa_map
                  + a.null_frac + misc,
              total_s);
  std::fflush(stdout);
  return total_s;
}

}  // namespace

int main(int argc, char** argv) {
  const int reps = (argc > 1) ? std::atoi(argv[1]) : 3;
  double t = 0.0;
  t += run_case("mcar30 C=6 R=5 n=40", 40, 5, 6, false, reps, 1);
  t += run_case("mcar30 C=6 R=5 n=100", 100, 5, 6, false, reps, 2);
  t += run_case("designed2 C=6 R=5 n=100", 100, 5, 6, true, reps, 3);
  t += run_case("mcar30 C=5 R=5 n=100", 100, 5, 5, false, reps, 4);
  t += run_case("designed2 C=6 R=5 n=300", 300, 5, 6, true, reps, 5);
  t += run_case("designed2 C=5 R=7 n=100", 100, 7, 5, true, 1, 6);
  t += run_case("flat1 designed2 C=6 R=5", 100, 5, 6, true, reps, 7, 1.0);
  t += run_case("flat1 mcar30 C=6 R=5", 100, 5, 6, false, reps, 8, 1.0);
  std::printf("grand total: %.2f s\n", t);
  return 0;
}
