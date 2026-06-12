#include "misskappa/estimate.hpp"
#include "misskappa/loss.hpp"
#include "prof.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using Clock = std::chrono::steady_clock;

namespace {

struct DgpSpec {
  std::string id;
  std::string label;
  std::string family;
  int C = 0;
  int R = 0;
  double agreement = 0.0;
  std::vector<double> category_prob;
};

struct Options {
  std::string variant = "unknown";
  std::string out = "results/replicates-unknown.csv";
  int reps = 5;
  int truth_n = 50000;
  int seed_base = 353500;
  std::vector<int> n_grid{40, 100};
  std::vector<std::string> dgps{"balanced6", "highagree6", "sparse6"};
  std::vector<std::string> mechanisms{"mcar30", "mar_anchor30", "designed_random2"};
  std::vector<std::string> weights{"identity", "quadratic"};
  bool progress = false;
};

void usage(const char* argv0) {
  std::cerr
      << "Usage: " << argv0 << " --variant LABEL --out PATH [options]\n\n"
      << "Options:\n"
      << "  --reps N              Replicates per cell. Default: 5.\n"
      << "  --truth-n N           Complete-data Monte Carlo truth size. Default: 50000.\n"
      << "  --seed-base N         Base seed. Default: 353500.\n"
      << "  --n-grid CSV          Sample sizes. Default: 40,100.\n"
      << "  --dgps CSV            DGP ids. Default: balanced6,highagree6,sparse6.\n"
      << "  --mechanisms CSV      Mechanism ids. Default: mcar30,mar_anchor30,designed_random2.\n"
      << "  --weights CSV         identity,quadratic. Default: identity,quadratic.\n"
      << "  --progress            Print one line per design cell.\n";
}

std::vector<std::string> split_csv(const std::string& text) {
  std::vector<std::string> out;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    item.erase(item.begin(), std::find_if(item.begin(), item.end(), [](unsigned char ch) {
      return !std::isspace(ch);
    }));
    item.erase(std::find_if(item.rbegin(), item.rend(), [](unsigned char ch) {
      return !std::isspace(ch);
    }).base(), item.end());
    if (!item.empty()) out.push_back(item);
  }
  return out;
}

std::vector<int> split_csv_int(const std::string& text) {
  std::vector<int> out;
  for (const std::string& item : split_csv(text)) out.push_back(std::atoi(item.c_str()));
  return out;
}

Options parse_args(int argc, char** argv) {
  Options opts;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << name << " needs a value.\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      std::exit(0);
    } else if (arg == "--variant") {
      opts.variant = value("--variant");
    } else if (arg == "--out") {
      opts.out = value("--out");
    } else if (arg == "--reps") {
      opts.reps = std::atoi(value("--reps").c_str());
    } else if (arg == "--truth-n") {
      opts.truth_n = std::atoi(value("--truth-n").c_str());
    } else if (arg == "--seed-base") {
      opts.seed_base = std::atoi(value("--seed-base").c_str());
    } else if (arg == "--n-grid") {
      opts.n_grid = split_csv_int(value("--n-grid"));
    } else if (arg == "--dgps") {
      opts.dgps = split_csv(value("--dgps"));
    } else if (arg == "--mechanisms") {
      opts.mechanisms = split_csv(value("--mechanisms"));
    } else if (arg == "--weights") {
      opts.weights = split_csv(value("--weights"));
    } else if (arg == "--progress") {
      opts.progress = true;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      std::exit(2);
    }
  }
  if (opts.reps < 1 || opts.truth_n < 1000 || opts.n_grid.empty() ||
      opts.dgps.empty() || opts.mechanisms.empty() || opts.weights.empty()) {
    std::cerr << "Invalid empty grid or non-positive replicate/truth size.\n";
    std::exit(2);
  }
  return opts;
}

DgpSpec dgp_by_id(const std::string& id) {
  if (id == "balanced6") {
    return DgpSpec{id, "balanced C6 R5", "latent_class", 6, 5, 0.70,
                   std::vector<double>(6, 1.0 / 6.0)};
  }
  if (id == "highagree6") {
    return DgpSpec{id, "high-agreement C6 R5", "latent_class", 6, 5, 0.88,
                   std::vector<double>(6, 1.0 / 6.0)};
  }
  if (id == "sparse6") {
    return DgpSpec{id, "skew-sparse C6 R5", "latent_class", 6, 5, 0.75,
                   {0.52, 0.19, 0.11, 0.08, 0.06, 0.04}};
  }
  if (id == "sparse5r7") {
    return DgpSpec{id, "skew-sparse C5 R7", "latent_class", 5, 7, 0.76,
                   {0.46, 0.22, 0.15, 0.10, 0.07}};
  }
  std::cerr << "Unknown DGP id: " << id << "\n";
  std::exit(2);
}

int draw_category(const std::vector<double>& p, std::mt19937& rng) {
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  const double u = uni(rng);
  double acc = 0.0;
  for (std::size_t k = 0; k < p.size(); ++k) {
    acc += p[k];
    if (u <= acc) return static_cast<int>(k);
  }
  return static_cast<int>(p.size() - 1);
}

int draw_category_excluding(const std::vector<double>& p, int excluded, std::mt19937& rng) {
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  const double total = std::max(1e-12, 1.0 - p[static_cast<std::size_t>(excluded)]);
  const double u = uni(rng) * total;
  double acc = 0.0;
  for (std::size_t k = 0; k < p.size(); ++k) {
    if (static_cast<int>(k) == excluded) continue;
    acc += p[k];
    if (u <= acc) return static_cast<int>(k);
  }
  return (excluded == 0) ? 1 : 0;
}

misskappa::IntMat simulate_complete(const DgpSpec& spec, int n, std::mt19937& rng) {
  misskappa::IntMat x(n, spec.R);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  for (int i = 0; i < n; ++i) {
    const int truth = draw_category(spec.category_prob, rng);
    for (int r = 0; r < spec.R; ++r) {
      x(i, r) = (uni(rng) < spec.agreement)
                    ? truth
                    : draw_category_excluding(spec.category_prob, truth, rng);
    }
  }
  return x;
}

void restore_min_observed(
    misskappa::IntMat& x, const misskappa::IntMat& full, int row, int min_obs,
    std::mt19937& rng) {
  std::uniform_int_distribution<int> pick_col(0, static_cast<int>(x.cols()) - 1);
  for (;;) {
    int obs = 0;
    std::vector<int> missing;
    for (int j = 0; j < x.cols(); ++j) {
      if (x(row, j) == misskappa::na_code) missing.push_back(j);
      else ++obs;
    }
    if (obs >= min_obs || missing.empty()) return;
    int j = missing[static_cast<std::size_t>(pick_col(rng) % missing.size())];
    x(row, j) = full(row, j);
  }
}

void apply_mcar(misskappa::IntMat& x, const misskappa::IntMat& full, double p, std::mt19937& rng) {
  std::bernoulli_distribution miss(p);
  for (int i = 0; i < x.rows(); ++i) {
    for (int j = 0; j < x.cols(); ++j) {
      if (miss(rng)) x(i, j) = misskappa::na_code;
    }
    restore_min_observed(x, full, i, 2, rng);
  }
}

void apply_mar_anchor30(misskappa::IntMat& x, const misskappa::IntMat& full, int C, std::mt19937& rng) {
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  for (int i = 0; i < x.rows(); ++i) {
    x(i, 0) = full(i, 0);
    const double anchor = (C > 1) ? static_cast<double>(full(i, 0)) / static_cast<double>(C - 1) : 0.0;
    const double p = std::min(0.75, std::max(0.03, 0.12 + 0.36 * anchor));
    for (int j = 1; j < x.cols(); ++j) {
      if (uni(rng) < p) x(i, j) = misskappa::na_code;
    }
    restore_min_observed(x, full, i, 2, rng);
  }
}

void comb_rec(
    int R, int k, int start, std::vector<int>& current,
    std::vector<std::vector<int>>& out) {
  if (static_cast<int>(current.size()) == k) {
    out.push_back(current);
    return;
  }
  const int need = k - static_cast<int>(current.size());
  for (int r = start; r <= R - need; ++r) {
    current.push_back(r);
    comb_rec(R, k, r + 1, current, out);
    current.pop_back();
  }
}

std::vector<std::vector<int>> combinations(int R, int k) {
  std::vector<std::vector<int>> out;
  std::vector<int> current;
  comb_rec(R, k, 0, current, out);
  return out;
}

void keep_only(misskappa::IntMat& x, int row, const std::vector<int>& keep) {
  for (int j = 0; j < x.cols(); ++j) x(row, j) = misskappa::na_code;
  for (int j : keep) x(row, j) = 0;  // placeholder restored by caller
}

void apply_designed(misskappa::IntMat& x, const misskappa::IntMat& full, int keep, std::mt19937& rng) {
  const int R = static_cast<int>(x.cols());
  std::vector<std::vector<int>> blocks = combinations(R, keep);
  std::uniform_int_distribution<std::size_t> pick(0, blocks.size() - 1);
  for (int i = 0; i < x.rows(); ++i) {
    const std::vector<int>& block =
        (i < static_cast<int>(blocks.size())) ? blocks[static_cast<std::size_t>(i)] : blocks[pick(rng)];
    for (int j = 0; j < R; ++j) x(i, j) = misskappa::na_code;
    for (int j : block) x(i, j) = full(i, j);
  }
}

std::string mechanism_family(const std::string& mechanism) {
  if (mechanism.find("mcar") == 0) return "mcar";
  if (mechanism.find("mar") == 0) return "mar";
  if (mechanism.find("designed") == 0) return "planned";
  return "unknown";
}

misskappa::IntMat apply_missing(
    const misskappa::IntMat& full, const DgpSpec& spec,
    const std::string& mechanism, std::mt19937& rng) {
  misskappa::IntMat x = full;
  if (mechanism == "complete") return x;
  if (mechanism == "mcar30") {
    apply_mcar(x, full, 0.30, rng);
    return x;
  }
  if (mechanism == "mar_anchor30") {
    apply_mar_anchor30(x, full, spec.C, rng);
    return x;
  }
  if (mechanism == "designed_random2") {
    apply_designed(x, full, 2, rng);
    return x;
  }
  if (mechanism == "designed_bib3") {
    apply_designed(x, full, 3, rng);
    return x;
  }
  std::cerr << "Unknown mechanism: " << mechanism << "\n";
  std::exit(2);
}

misskappa::RealMat make_weights(const std::string& weight, int C) {
  if (weight == "identity" || weight == "nominal") {
    auto W = misskappa::loss::identity_weights(C);
    if (!W) std::exit(2);
    return *W;
  }
  if (weight == "quadratic") {
    misskappa::RealVec values(C);
    for (int k = 0; k < C; ++k) values(k) = static_cast<double>(k);
    auto W = misskappa::loss::quadratic_weights(C, values);
    if (!W) std::exit(2);
    return *W;
  }
  std::cerr << "Unknown weight: " << weight << "\n";
  std::exit(2);
}

std::string weight_label(const std::string& weight) {
  return (weight == "identity" || weight == "nominal") ? "nominal" : weight;
}

std::string error_name(misskappa::Error e) {
  switch (e) {
    case misskappa::Error::invalid_argument: return "invalid_argument";
    case misskappa::Error::dimension_mismatch: return "dimension_mismatch";
    case misskappa::Error::singular_weight: return "singular_weight";
    case misskappa::Error::numerical_error: return "numerical_error";
    case misskappa::Error::not_identified: return "not_identified";
    case misskappa::Error::not_supported: return "not_supported";
    case misskappa::Error::not_converged: return "not_converged";
  }
  return "unknown";
}

int complete_rows(const misskappa::IntMat& x) {
  int out = 0;
  for (int i = 0; i < x.rows(); ++i) {
    bool complete = true;
    for (int j = 0; j < x.cols(); ++j) complete = complete && x(i, j) != misskappa::na_code;
    if (complete) ++out;
  }
  return out;
}

int empty_rows(const misskappa::IntMat& x) {
  int out = 0;
  for (int i = 0; i < x.rows(); ++i) {
    bool empty = true;
    for (int j = 0; j < x.cols(); ++j) empty = empty && x(i, j) == misskappa::na_code;
    if (empty) ++out;
  }
  return out;
}

double observed_fraction(const misskappa::IntMat& x) {
  int obs = 0;
  for (int i = 0; i < x.rows(); ++i)
    for (int j = 0; j < x.cols(); ++j)
      if (x(i, j) != misskappa::na_code) ++obs;
  return static_cast<double>(obs) / static_cast<double>(x.rows() * x.cols());
}

int min_pair_count(const misskappa::IntMat& x) {
  int ans = x.rows();
  for (int a = 0; a < x.cols() - 1; ++a) {
    for (int b = a + 1; b < x.cols(); ++b) {
      int count = 0;
      for (int i = 0; i < x.rows(); ++i) {
        if (x(i, a) != misskappa::na_code && x(i, b) != misskappa::na_code) ++count;
      }
      ans = std::min(ans, count);
    }
  }
  return ans;
}

int observed_patterns(const misskappa::IntMat& x) {
  std::set<std::string> patterns;
  for (int i = 0; i < x.rows(); ++i) {
    std::string key;
    for (int j = 0; j < x.cols(); ++j) {
      key += std::to_string(x(i, j));
      key += ':';
    }
    patterns.insert(std::move(key));
  }
  return static_cast<int>(patterns.size());
}

void write_num(std::ostream& os, double x) {
  if (std::isfinite(x)) os << std::setprecision(17) << x;
  else os << "NA";
}

void write_header(std::ostream& os) {
  os << "variant,method,estimator,dgp,dgp_label,dgp_family,C,R,mechanism,"
        "mechanism_family,n,rep,seed,weight,weight_label,coefficient,truth,"
        "estimate,se,null_frac,n_eff,error,elapsed_ms,em_iters,support,active,"
        "groups,observed_fraction,subjects_used,empty_rows,complete_rows,"
        "min_pair_count,observed_patterns\n";
}

void write_row(
    std::ostream& os, const Options& opts, const DgpSpec& spec,
    const std::string& mechanism, int n, int rep, int seed,
    const std::string& weight, const std::string& coefficient, double truth,
    double estimate, double se, double null_frac, int n_eff,
    const std::string& error, double elapsed_ms, const bench_prof::Accum& acc,
    double obs_frac, int subjects_used, int empty, int complete, int min_pair,
    int obs_patterns) {
  os << opts.variant << ",cat_fiml_" << opts.variant << ",cat_fiml,"
     << spec.id << "," << spec.label << "," << spec.family << ","
     << spec.C << "," << spec.R << "," << mechanism << ","
     << mechanism_family(mechanism) << "," << n << "," << rep << ","
     << seed << "," << weight << "," << weight_label(weight) << ","
     << coefficient << ",";
  write_num(os, truth);
  os << ",";
  write_num(os, estimate);
  os << ",";
  write_num(os, se);
  os << ",";
  write_num(os, null_frac);
  os << "," << n_eff << "," << error << ",";
  write_num(os, elapsed_ms);
  os << "," << acc.em_iters << "," << acc.n_final << "," << acc.n_active
     << "," << acc.n_groups << ",";
  write_num(os, obs_frac);
  os << "," << subjects_used << "," << empty << "," << complete << ","
     << min_pair << "," << obs_patterns << "\n";
}

std::vector<double> truth_for(const DgpSpec& spec, const std::string& weight, int truth_n, int seed) {
  std::mt19937 rng(static_cast<std::uint32_t>(seed));
  misskappa::IntMat x = simulate_complete(spec, truth_n, rng);
  misskappa::RealMat W = make_weights(weight, spec.C);
  auto fit = misskappa::estimate_available(x, W);
  if (!fit) {
    std::cerr << "Truth fit failed for " << spec.id << " " << weight << ": "
              << error_name(fit.error()) << "\n";
    std::exit(2);
  }
  return {fit->estimates(0), fit->estimates(1), fit->estimates(2)};
}

}  // namespace

int main(int argc, char** argv) {
  const Options opts = parse_args(argc, argv);
  std::ofstream out(opts.out);
  if (!out) {
    std::cerr << "Failed to open output file: " << opts.out << "\n";
    return 2;
  }
  write_header(out);

  const std::vector<std::string> coef_names{"Conger", "Fleiss", "Brennan_Prediger"};
  misskappa::EmOptions em_opts;
  em_opts.tol = 1e-7;
  em_opts.max_iter = 12000;
  em_opts.prune_tol = 1e-10;
  em_opts.start_alpha = 0.1;
  em_opts.info_rcond = 1e-4;
  em_opts.flatten = 0.0;

  for (const std::string& dgp_id : opts.dgps) {
    const DgpSpec spec = dgp_by_id(dgp_id);
    for (const std::string& weight : opts.weights) {
      const int truth_seed = opts.seed_base + 700000 + static_cast<int>(31 * spec.C + 17 * spec.R)
                             + (weight == "quadratic" ? 1009 : 0);
      const std::vector<double> truth = truth_for(spec, weight, opts.truth_n, truth_seed);
      const misskappa::RealMat W = make_weights(weight, spec.C);

      for (const std::string& mechanism : opts.mechanisms) {
        for (const int n : opts.n_grid) {
          if (opts.progress) {
            std::cerr << opts.variant << ": " << spec.id << " " << mechanism
                      << " n=" << n << " weight=" << weight
                      << " reps=" << opts.reps << "\n";
          }
          for (int rep = 1; rep <= opts.reps; ++rep) {
            const int seed = opts.seed_base + 1000000 * (rep + 13 * n)
                             + 1000 * static_cast<int>(mechanism.size())
                             + 101 * spec.C + spec.R;
            std::mt19937 rng(static_cast<std::uint32_t>(seed));
            const misskappa::IntMat full = simulate_complete(spec, n, rng);
            misskappa::IntMat x = apply_missing(full, spec, mechanism, rng);
            const double obs_frac = observed_fraction(x);
            const int empty = empty_rows(x);
            const int complete = complete_rows(x);
            const int min_pair = min_pair_count(x);
            const int obs_patterns = observed_patterns(x);

            bench_prof::acc().reset();
            const auto t0 = Clock::now();
            auto fit = misskappa::estimate_fiml(x, W, em_opts);
            const auto t1 = Clock::now();
            const double elapsed_ms =
                1000.0 * std::chrono::duration<double>(t1 - t0).count();
            const bench_prof::Accum acc = bench_prof::acc();

            if (!fit) {
              const std::string err = error_name(fit.error());
              for (std::size_t k = 0; k < coef_names.size(); ++k) {
                write_row(out, opts, spec, mechanism, n, rep, seed, weight,
                          coef_names[k], truth[k], NAN, NAN, NAN, n, err,
                          elapsed_ms, acc, obs_frac, n, empty, complete,
                          min_pair, obs_patterns);
              }
              continue;
            }

            for (std::size_t k = 0; k < coef_names.size(); ++k) {
              double se = NAN;
              if (fit->vcov.rows() > static_cast<Eigen::Index>(k) &&
                  fit->vcov.cols() > static_cast<Eigen::Index>(k)) {
                const double var = fit->vcov(static_cast<Eigen::Index>(k),
                                             static_cast<Eigen::Index>(k));
                if (std::isfinite(var) && var >= 0.0) se = std::sqrt(var);
              }
              const double nf =
                  (fit->null_frac.size() > static_cast<Eigen::Index>(k))
                      ? fit->null_frac(static_cast<Eigen::Index>(k))
                      : NAN;
              write_row(out, opts, spec, mechanism, n, rep, seed, weight,
                        coef_names[k], truth[k],
                        fit->estimates(static_cast<Eigen::Index>(k)), se, nf,
                        n, "", elapsed_ms, acc, obs_frac, n, empty, complete,
                        min_pair, obs_patterns);
            }
          }
        }
      }
    }
  }

  return 0;
}
