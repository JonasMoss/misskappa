#!/usr/bin/env Rscript
#
# Experiment 29: finite-sample agreement under accidental and planned
# missing ratings. The grid extends the coverage runner with small-n cells,
# planned-missingness designs, listwise deletion, Moss/Perreault-Leigh
# high-agreement cells, ordinal normal-copula data, and Dawid-Skene-style
# severity/confusion-matrix raters.

suppressPackageStartupMessages({
  library(misskappa)
})

usage <- function(status = 0L) {
  cat(paste0(
    "Usage: Rscript run_experiment.R [mode/options]\n\n",
    "Modes:\n",
    "  --smoke            Mechanical check: n=20, reps=5, quadratic methods.\n",
    "  --screen           Broad screen: n=10,20,40,100,300, reps=500.\n",
    "  --targeted         Paper-grade targeted defaults: reps=2000.\n\n",
    "Options:\n",
    "  --reps N            Override replications per cell.\n",
    "  --n-grid CSV        Override sample sizes.\n",
    "  --truth-n N         Complete-data Monte Carlo truth size.\n",
    "  --dgps CSV          DGP ids.\n",
    "  --mechanisms CSV    Missingness mechanism ids.\n",
    "  --methods CSV       Method ids.\n",
    "  --seed-base N       Base seed. Default: 292900.\n",
    "  --out-dir PATH      Output directory. Default: script-local results/.\n",
    "  --shard-index N     1-based shard index. Default: 1.\n",
    "  --shard-count N     Number of replicate shards. Default: 1.\n",
    "  --resume            Reuse existing per-cell checkpoints in --out-dir.\n",
    "  --progress          Print one line per design cell.\n",
    "  --help, -h          Show this help and exit.\n\n",
    "Examples:\n",
    "  Rscript run_experiment.R --smoke --progress\n",
    "  Rscript run_experiment.R --screen --shard-index 1 --shard-count 8 --progress\n",
    "  Rscript run_experiment.R --screen --reps 3000 --resume --progress\n"
  ))
  quit(save = "no", status = status)
}

script_dir <- local({
  cmd <- commandArgs(FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0L) getwd()
  else dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
})

parse_csv_chr <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_csv_int <- function(x, arg) {
  out <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (any(is.na(out))) stop(arg, " must be a comma-separated integer list.", call. = FALSE)
  unique(out)
}

all_methods <- c(
  "ipw_nominal", "cat_fiml_nominal", "listwise_nominal",
  "ipw_absolute", "cat_fiml_absolute", "listwise_absolute",
  "ipw_quadratic", "cat_fiml_quadratic", "pairwise_quadratic",
  "nt_fiml_quadratic", "listwise_quadratic"
)

quadratic_methods <- c(
  "ipw_quadratic", "cat_fiml_quadratic", "pairwise_quadratic",
  "nt_fiml_quadratic", "listwise_quadratic"
)

mode_defaults <- function(mode) {
  if (mode == "smoke") {
    return(list(
      mode = mode,
      reps = 5L,
      n_grid = 20L,
      truth_n = 3000L,
      dgps = c("balanced4", "skill5x6_sparse", "moss5x5_k80",
               "copula5x6_high", "ds5x6_severity"),
      mechanisms = c("complete", "mcar30", "mar_anchor30",
                     "designed_bib3", "designed_two_phase"),
      methods = quadratic_methods
    ))
  }
  if (mode == "screen") {
    return(list(
      mode = mode,
      reps = 500L,
      n_grid = c(10L, 20L, 40L, 100L, 300L),
      truth_n = 150000L,
      dgps = c("balanced4", "biased4", "sparse5", "skill5x6_sparse",
               "moss5x2_k50", "moss5x2_k80", "moss5x2_k90",
               "moss5x5_k50", "moss5x5_k80", "moss5x5_k90",
               "copula5x6_mid", "copula5x6_high", "copula5x6_skew",
               "ds5x6_symmetric", "ds5x6_severity"),
      mechanisms = c("complete", "mcar15", "mcar30", "mcar_rater30",
                     "mcar_subject30", "mar_anchor15", "mar_anchor30",
                     "mar_anchor_nonlinear30", "mar_shifted30",
                     "designed_bib3", "designed_forms4",
                     "designed_anchor3", "designed_two_phase"),
      methods = all_methods
    ))
  }
  if (mode == "targeted") {
    return(list(
      mode = mode,
      reps = 2000L,
      n_grid = c(10L, 20L, 40L, 100L),
      truth_n = 300000L,
      dgps = c("skill5x6_sparse", "moss5x5_k80", "moss5x5_k90",
               "copula5x6_high", "ds5x6_severity"),
      mechanisms = c("complete", "mcar30", "mar_anchor30",
                     "designed_bib3", "designed_anchor3",
                     "designed_two_phase"),
      methods = all_methods
    ))
  }
  stop("Unknown mode: ", mode, call. = FALSE)
}

parse_args <- function(argv) {
  mode <- "smoke"
  if ("--smoke" %in% argv) mode <- "smoke"
  if ("--screen" %in% argv) mode <- "screen"
  if ("--targeted" %in% argv) mode <- "targeted"
  opts <- mode_defaults(mode)
  opts$seed_base <- 292900L
  opts$out_dir <- file.path(script_dir, "results")
  opts$progress <- FALSE
  opts$resume <- FALSE
  opts$shard_index <- 1L
  opts$shard_count <- 1L

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg %in% c("--help", "-h")) usage(0L)
    if (arg %in% c("--smoke", "--screen", "--targeted")) {
      i <- i + 1L
      next
    }
    if (arg == "--progress") {
      opts$progress <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "--resume") {
      opts$resume <- TRUE
      i <- i + 1L
      next
    }
    needs_value <- c("--reps", "--n-grid", "--truth-n", "--dgps",
                     "--mechanisms", "--methods", "--seed-base", "--out-dir",
                     "--shard-index", "--shard-count")
    if (arg %in% needs_value) {
      if (i == length(argv)) stop(arg, " needs a value.", call. = FALSE)
      val <- argv[[i + 1L]]
      if (arg == "--reps") opts$reps <- as.integer(val)
      else if (arg == "--n-grid") opts$n_grid <- parse_csv_int(val, arg)
      else if (arg == "--truth-n") opts$truth_n <- as.integer(val)
      else if (arg == "--dgps") opts$dgps <- parse_csv_chr(val)
      else if (arg == "--mechanisms") opts$mechanisms <- parse_csv_chr(val)
      else if (arg == "--methods") opts$methods <- parse_csv_chr(val)
      else if (arg == "--seed-base") opts$seed_base <- as.integer(val)
      else if (arg == "--out-dir") opts$out_dir <- val
      else if (arg == "--shard-index") opts$shard_index <- as.integer(val)
      else if (arg == "--shard-count") opts$shard_count <- as.integer(val)
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", arg, call. = FALSE)
  }

  if (is.na(opts$reps) || opts$reps < 1L) stop("--reps must be >= 1.", call. = FALSE)
  if (any(is.na(opts$n_grid)) || any(opts$n_grid < 5L)) {
    stop("--n-grid entries must be integers >= 5.", call. = FALSE)
  }
  if (is.na(opts$truth_n) || opts$truth_n < 1000L) {
    stop("--truth-n must be an integer >= 1000.", call. = FALSE)
  }
  if (is.na(opts$seed_base)) stop("--seed-base must be an integer.", call. = FALSE)
  if (is.na(opts$shard_count) || opts$shard_count < 1L) {
    stop("--shard-count must be >= 1.", call. = FALSE)
  }
  if (is.na(opts$shard_index) || opts$shard_index < 1L ||
      opts$shard_index > opts$shard_count) {
    stop("--shard-index must be between 1 and --shard-count.", call. = FALSE)
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

# ---- DGPs -------------------------------------------------------------------

centered_guess <- function(C, centers, spread = 0.9) {
  vals <- seq_len(C)
  G <- vapply(centers, function(ctr) {
    p <- exp(-0.5 * ((vals - ctr) / spread)^2)
    p / sum(p)
  }, numeric(C))
  t(G)
}

skill_spec <- function(id, label, C, R, p, skill, guess, truth_skill_mult) {
  list(
    id = id,
    label = label,
    family = "skill_guess",
    type = "skill",
    C = as.integer(C),
    R = as.integer(R),
    p = p,
    skill = skill,
    guess = guess,
    truth_skill_mult = truth_skill_mult
  )
}

moss_spec <- function(R, target) {
  list(
    id = sprintf("moss5x%d_k%d", R, round(100 * target)),
    label = sprintf("Moss/Perreault-Leigh five-category R=%d target %.2f", R, target),
    family = "moss_pl",
    type = "moss",
    C = 5L,
    R = as.integer(R),
    p = rep(0.2, 5L),
    target_quadratic = target,
    knowledge_prob = sqrt(target)
  )
}

copula_spec <- function(id, label, corr, p) {
  list(
    id = id,
    label = label,
    family = "ordinal_copula",
    type = "copula",
    C = 5L,
    R = 6L,
    p = p,
    latent_corr = corr
  )
}

ds_spec <- function(id, label, severity) {
  list(
    id = id,
    label = label,
    family = "dawid_skene_style",
    type = "ds",
    C = 5L,
    R = 6L,
    p = c(0.38, 0.27, 0.18, 0.11, 0.06),
    accuracy = c(0.86, 0.82, 0.78, 0.74, 0.70, 0.66),
    severity = severity
  )
}

dgp_defs <- local({
  out <- list(
    balanced4 = skill_spec(
      "balanced4",
      "Balanced four-category, four-rater skill model",
      4L, 4L, rep(0.25, 4L), rep(0.72, 4L),
      matrix(rep(0.25, 16L), nrow = 4L, byrow = TRUE),
      rep(1, 4L)
    ),
    biased4 = skill_spec(
      "biased4",
      "Biased four-category, non-exchangeable raters",
      4L, 4L, c(0.45, 0.30, 0.18, 0.07),
      c(0.88, 0.78, 0.66, 0.58),
      centered_guess(4L, c(1.2, 1.8, 3.0, 3.7), spread = 0.75),
      c(1.05, 0.98, 0.86, 0.70)
    ),
    sparse5 = skill_spec(
      "sparse5",
      "Sparse five-category, five-rater high-agreement model",
      5L, 5L, c(0.55, 0.22, 0.12, 0.07, 0.04),
      c(0.90, 0.86, 0.82, 0.78, 0.74),
      centered_guess(5L, c(1.1, 1.7, 2.6, 3.7, 4.5), spread = 0.8),
      c(1.05, 0.98, 0.88, 0.74, 0.58)
    ),
    skill5x6_sparse = skill_spec(
      "skill5x6_sparse",
      "Sparse five-category, six-rater Cat-FIML stress model",
      5L, 6L, c(0.50, 0.22, 0.13, 0.09, 0.06),
      c(0.88, 0.84, 0.80, 0.76, 0.72, 0.68),
      centered_guess(5L, c(1.1, 1.5, 2.2, 3.0, 4.0, 4.6), spread = 0.8),
      c(1.06, 0.99, 0.89, 0.76, 0.60)
    ),
    copula5x6_mid = copula_spec(
      "copula5x6_mid",
      "Five-category six-rater ordinal normal-copula, mid agreement",
      0.45, rep(0.2, 5L)
    ),
    copula5x6_high = copula_spec(
      "copula5x6_high",
      "Five-category six-rater ordinal normal-copula, high agreement",
      0.72, rep(0.2, 5L)
    ),
    copula5x6_skew = copula_spec(
      "copula5x6_skew",
      "Five-category six-rater ordinal normal-copula, skewed marginals",
      0.55, c(0.44, 0.25, 0.16, 0.10, 0.05)
    ),
    ds5x6_symmetric = ds_spec(
      "ds5x6_symmetric",
      "Five-category six-rater Dawid-Skene-style symmetric errors",
      rep(0, 6L)
    ),
    ds5x6_severity = ds_spec(
      "ds5x6_severity",
      "Five-category six-rater Dawid-Skene-style severity shifts",
      c(-0.65, -0.35, -0.10, 0.15, 0.40, 0.70)
    )
  )
  for (R in c(2L, 5L, 6L)) {
    for (target in c(0.50, 0.80, 0.90)) {
      spec <- moss_spec(R, target)
      out[[spec$id]] <- spec
    }
  }
  out
})

simulate_skill <- function(n, spec) {
  truth <- sample.int(spec$C, n, replace = TRUE, prob = spec$p)
  X <- matrix(NA_integer_, nrow = n, ncol = spec$R)
  for (j in seq_len(spec$R)) {
    p_correct <- pmin(pmax(spec$skill[j] * spec$truth_skill_mult[truth], 0), 1)
    correct <- stats::rbinom(n, 1L, p_correct) == 1L
    guessed <- sample.int(spec$C, n, replace = TRUE, prob = spec$guess[j, ])
    X[, j] <- ifelse(correct, truth, guessed)
  }
  X
}

simulate_moss <- function(n, spec) {
  truth <- sample.int(spec$C, n, replace = TRUE)
  X <- matrix(NA_integer_, nrow = n, ncol = spec$R)
  for (j in seq_len(spec$R)) {
    knows <- stats::runif(n) < spec$knowledge_prob
    guesses <- sample.int(spec$C, n, replace = TRUE)
    X[, j] <- ifelse(knows, truth, guesses)
  }
  X
}

simulate_copula <- function(n, spec) {
  Sigma <- matrix(spec$latent_corr, spec$R, spec$R)
  diag(Sigma) <- 1
  Z <- matrix(stats::rnorm(n * spec$R), nrow = n, ncol = spec$R) %*% chol(Sigma)
  cuts <- c(-Inf, stats::qnorm(cumsum(spec$p)[-spec$C]), Inf)
  X <- apply(Z, 2L, function(z) as.integer(cut(z, breaks = cuts, labels = FALSE)))
  matrix(X, nrow = n, ncol = spec$R)
}

ds_row_prob <- function(C, truth, accuracy, severity) {
  cats <- seq_len(C)
  center <- truth + severity
  error <- exp(-0.5 * ((cats - center) / 0.85)^2)
  error[truth] <- 0
  if (sum(error) <= 0) error <- rep(1, C)
  error <- error / sum(error)
  p <- (1 - accuracy) * error
  p[truth] <- p[truth] + accuracy
  p / sum(p)
}

simulate_ds <- function(n, spec) {
  truth <- sample.int(spec$C, n, replace = TRUE, prob = spec$p)
  X <- matrix(NA_integer_, nrow = n, ncol = spec$R)
  for (j in seq_len(spec$R)) {
    for (i in seq_len(n)) {
      p <- ds_row_prob(spec$C, truth[i], spec$accuracy[j], spec$severity[j])
      X[i, j] <- sample.int(spec$C, 1L, prob = p)
    }
  }
  X
}

simulate_ratings <- function(n, spec) {
  X <- switch(
    spec$type,
    skill = simulate_skill(n, spec),
    moss = simulate_moss(n, spec),
    copula = simulate_copula(n, spec),
    ds = simulate_ds(n, spec),
    stop("Unknown DGP type: ", spec$type, call. = FALSE)
  )
  colnames(X) <- paste0("r", seq_len(spec$R))
  X
}

category_z <- function(x, C) {
  z <- (as.numeric(x) - mean(seq_len(C))) / stats::sd(seq_len(C))
  z[!is.finite(z)] <- 0
  z
}

tune_intercept <- function(eta, target) {
  f <- function(a) mean(stats::plogis(a + eta)) - target
  stats::uniroot(f, c(-20, 20))$root
}

drop_with_prob <- function(X, P, keep_first = FALSE) {
  out <- X
  M <- matrix(stats::runif(length(X)), nrow = nrow(X), ncol = ncol(X)) < P
  if (keep_first && ncol(out) >= 1L) M[, 1L] <- FALSE
  out[M] <- NA_integer_
  out
}

apply_mcar <- function(X, rate) {
  drop_with_prob(X, matrix(rate, nrow(X), ncol(X)))
}

apply_mcar_rater <- function(X, target = 0.30) {
  R <- ncol(X)
  rates <- seq(max(0.04, target - 0.22), min(0.70, target + 0.22), length.out = R)
  rates <- pmin(pmax(rates - mean(rates) + target, 0.02), 0.80)
  drop_with_prob(X, matrix(rates, nrow(X), R, byrow = TRUE))
}

apply_mcar_subject <- function(X, target = 0.30) {
  n <- nrow(X)
  p_i <- stats::rgamma(n, shape = 2.5, rate = 2.5)
  p_i <- pmin(0.85, p_i * target / mean(p_i))
  drop_with_prob(X, matrix(p_i, n, ncol(X)))
}

apply_anchor_mar <- function(X, C, target = 0.30, slope = 1.15,
                             nonlinear = FALSE, shifted = FALSE) {
  n <- nrow(X); R <- ncol(X)
  anchor <- category_z(X[, 1L], C)
  eta_base <- if (nonlinear) {
    centered <- anchor^2 - mean(anchor^2)
    1.35 * centered
  } else {
    slope * anchor
  }
  P <- matrix(0, n, R)
  for (j in seq_len(R)) {
    if (j == 1L) next
    shift <- if (shifted) (j - (R + 1) / 2) * 0.35 else 0
    eta <- eta_base + shift
    a <- tune_intercept(eta, target)
    P[, j] <- stats::plogis(a + eta)
  }
  drop_with_prob(X, P, keep_first = TRUE)
}

apply_sequential_mar <- function(X, C, target = 0.30) {
  out <- X
  R <- ncol(out)
  for (j in 2:R) {
    prev <- out[, j - 1L]
    prev_missing <- is.na(prev)
    z <- category_z(prev, C)
    eta <- 1.0 * z + ifelse(prev_missing, 1.25, 0)
    a <- tune_intercept(eta, target)
    p_miss <- stats::plogis(a + eta)
    out[stats::runif(nrow(out)) < p_miss, j] <- NA_integer_
  }
  out
}

cycle_patterns <- function(n, patterns) {
  patterns[((seq_len(n) - 1L) %% length(patterns)) + 1L]
}

apply_patterns <- function(X, patterns) {
  out <- X
  assigned <- cycle_patterns(nrow(X), patterns)
  for (i in seq_len(nrow(X))) {
    drop <- setdiff(seq_len(ncol(X)), assigned[[i]])
    if (length(drop)) out[i, drop] <- NA_integer_
  }
  out
}

apply_designed_bib <- function(X, k = 3L) {
  R <- ncol(X)
  k <- min(k, R)
  if (k >= R) return(X)
  patterns <- combn(seq_len(R), k, simplify = FALSE)
  apply_patterns(X, patterns)
}

apply_designed_forms4 <- function(X) {
  R <- ncol(X)
  if (R <= 4L) return(X)
  if (R >= 6L) {
    patterns <- list(c(1L, 2L, 3L, 4L), c(1L, 2L, 5L, 6L), c(3L, 4L, 5L, 6L))
  } else {
    patterns <- combn(seq_len(R), 4L, simplify = FALSE)
  }
  apply_patterns(X, patterns)
}

apply_designed_anchor3 <- function(X) {
  R <- ncol(X)
  if (R <= 3L) return(X)
  pairs <- combn(2:R, min(2L, R - 1L), simplify = FALSE)
  patterns <- lapply(pairs, function(p) c(1L, p))
  apply_patterns(X, patterns)
}

apply_designed_two_phase <- function(X) {
  out <- apply_designed_anchor3(X)
  full <- seq(5L, nrow(X), by = 5L)
  if (length(full)) out[full, ] <- X[full, ]
  out
}

apply_missing <- function(X, spec, mechanism) {
  if (mechanism == "complete") return(X)
  if (mechanism == "mcar15") return(apply_mcar(X, 0.15))
  if (mechanism == "mcar30") return(apply_mcar(X, 0.30))
  if (mechanism == "mcar_rater30") return(apply_mcar_rater(X, 0.30))
  if (mechanism == "mcar_subject30") return(apply_mcar_subject(X, 0.30))
  if (mechanism == "mar_anchor15") return(apply_anchor_mar(X, spec$C, 0.15))
  if (mechanism == "mar_anchor30") return(apply_anchor_mar(X, spec$C, 0.30))
  if (mechanism == "mar_anchor_nonlinear30") {
    return(apply_anchor_mar(X, spec$C, 0.30, nonlinear = TRUE))
  }
  if (mechanism == "mar_shifted30") {
    return(apply_anchor_mar(X, spec$C, 0.30, shifted = TRUE))
  }
  if (mechanism == "mar_sequential30") return(apply_sequential_mar(X, spec$C, 0.30))
  if (mechanism == "designed_bib3") return(apply_designed_bib(X, 3L))
  if (mechanism == "designed_forms4") return(apply_designed_forms4(X))
  if (mechanism == "designed_anchor3") return(apply_designed_anchor3(X))
  if (mechanism == "designed_two_phase") return(apply_designed_two_phase(X))
  stop("Unknown mechanism: ", mechanism, call. = FALSE)
}

# ---- estimators and intervals ----------------------------------------------

method_defs <- data.frame(
  method = all_methods,
  estimator = c("ipw", "cat_fiml", "listwise",
                "ipw", "cat_fiml", "listwise",
                "ipw", "cat_fiml", "pairwise", "nt_fiml", "listwise"),
  base_estimator = c("ipw", "cat_fiml", "ipw",
                     "ipw", "cat_fiml", "ipw",
                     "ipw", "cat_fiml", "pairwise", "nt_fiml", "ipw"),
  weight = c("nominal", "nominal", "nominal",
             "linear", "linear", "linear",
             "quadratic", "quadratic", "quadratic", "quadratic", "quadratic"),
  weight_label = c("nominal", "nominal", "nominal",
                   "absolute", "absolute", "absolute",
                   "quadratic", "quadratic", "quadratic", "quadratic", "quadratic"),
  stringsAsFactors = FALSE
)

valid_dgps <- names(dgp_defs)
valid_mechanisms <- c("complete", "mcar15", "mcar30", "mcar_rater30",
                      "mcar_subject30", "mar_anchor15", "mar_anchor30",
                      "mar_anchor_nonlinear30", "mar_shifted30",
                      "mar_sequential30", "designed_bib3", "designed_forms4",
                      "designed_anchor3", "designed_two_phase")
valid_methods <- method_defs$method

if (!all(opts$dgps %in% valid_dgps)) {
  stop("--dgps must contain only: ", paste(valid_dgps, collapse = ","), call. = FALSE)
}
if (!all(opts$mechanisms %in% valid_mechanisms)) {
  stop("--mechanisms must contain only: ", paste(valid_mechanisms, collapse = ","), call. = FALSE)
}
if (!all(opts$methods %in% valid_methods)) {
  stop("--methods must contain only: ", paste(valid_methods, collapse = ","), call. = FALSE)
}
method_defs <- method_defs[match(opts$methods, method_defs$method), , drop = FALSE]

coef_names <- c("Conger", "Fleiss")
alpha <- 0.05

se_vector <- function(fit, coefs) {
  V <- stats::vcov(fit)
  out <- rep(NA_real_, length(coefs)); names(out) <- coefs
  have <- intersect(coefs, rownames(V))
  vals <- diag(V)[have]
  vals[!is.finite(vals) | vals < 0] <- NA_real_
  out[have] <- sqrt(vals)
  out
}

clip_unit <- function(x) pmin(pmax(x, -1 + 1e-10), 1 - 1e-10)

make_intervals <- function(est, se, n_eff) {
  names <- c("wald_z", "wald_t", "fisher_t", "asin_t")
  out <- matrix(NA_real_, nrow = length(names), ncol = 2L,
                dimnames = list(names, c("lower", "upper")))
  if (!is.finite(est) || !is.finite(se) || se < 0 || !is.finite(n_eff) || n_eff < 2) {
    return(out)
  }
  z <- stats::qnorm(1 - alpha / 2)
  tcrit <- stats::qt(1 - alpha / 2, df = max(1, n_eff - 1L))
  out["wald_z", ] <- c(est - z * se, est + z * se)
  out["wald_t", ] <- c(est - tcrit * se, est + tcrit * se)

  est_c <- clip_unit(est)
  se_f <- se / (1 - est_c^2)
  g_f <- atanh(est_c)
  out["fisher_t", ] <- tanh(c(g_f - tcrit * se_f, g_f + tcrit * se_f))

  se_a <- se / sqrt(1 - est_c^2)
  g_a <- asin(est_c)
  out["asin_t", ] <- sin(c(g_a - tcrit * se_a, g_a + tcrit * se_a))
  out
}

em_options_for <- function(estimator) {
  if (estimator == "cat_fiml") {
    return(list(tol = 1e-7, max_iter = 12000L, prune_tol = 1e-10,
                start_alpha = 0.1, info_rcond = 1e-4))
  }
  if (estimator == "nt_fiml") {
    return(list(tol = 1e-8, max_iter = 12000L, fd_h = 1e-5))
  }
  list()
}

score_values <- function(X) {
  cats <- sort(unique(as.integer(X[!is.na(X)])))
  if (length(cats) < 2L) stop("Fewer than two observed categories.", call. = FALSE)
  cats - 1L
}

fit_kappa <- function(X, spec, estimator, weight) {
  misskappa::kappa(
    X,
    estimator = estimator,
    weight = weight,
    values = score_values(X),
    em_options = em_options_for(estimator)
  )
}

fit_one <- function(X, spec, method_row) {
  start <- proc.time()[["elapsed"]]
  out <- tryCatch({
    X_fit <- X
    if (method_row$estimator == "listwise") {
      X_fit <- X[stats::complete.cases(X), , drop = FALSE]
      if (nrow(X_fit) < 2L) stop("Listwise deletion left fewer than two complete rows.", call. = FALSE)
    }
    fit <- fit_kappa(X_fit, spec, method_row$base_estimator, method_row$weight)
    est <- stats::coef(fit)
    se <- se_vector(fit, coef_names)
    n_eff <- nrow(X_fit)

    rows <- lapply(coef_names, function(coef) {
      ci <- make_intervals(est[[coef]], se[[coef]], n_eff)
      data.frame(
        coefficient = coef,
        estimate = unname(est[[coef]]),
        se = unname(se[[coef]]),
        n_eff = n_eff,
        wald_z_lower = ci["wald_z", "lower"],
        wald_z_upper = ci["wald_z", "upper"],
        wald_t_lower = ci["wald_t", "lower"],
        wald_t_upper = ci["wald_t", "upper"],
        fisher_t_lower = ci["fisher_t", "lower"],
        fisher_t_upper = ci["fisher_t", "upper"],
        asin_t_lower = ci["asin_t", "lower"],
        asin_t_upper = ci["asin_t", "upper"],
        error = "",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }, error = function(e) {
    data.frame(
      coefficient = coef_names,
      estimate = NA_real_,
      se = NA_real_,
      n_eff = NA_integer_,
      wald_z_lower = NA_real_,
      wald_z_upper = NA_real_,
      wald_t_lower = NA_real_,
      wald_t_upper = NA_real_,
      fisher_t_lower = NA_real_,
      fisher_t_upper = NA_real_,
      asin_t_lower = NA_real_,
      asin_t_upper = NA_real_,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  out$elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start)
  out
}

min_pair_count <- function(X) {
  R <- ncol(X)
  if (R < 2L) return(0L)
  counts <- integer(R * (R - 1L) / 2L)
  idx <- 0L
  for (j in seq_len(R - 1L)) {
    for (k in (j + 1L):R) {
      idx <- idx + 1L
      counts[idx] <- sum(!is.na(X[, j]) & !is.na(X[, k]))
    }
  }
  min(counts)
}

observed_pattern_count <- function(X) {
  keys <- apply(X, 1L, function(row) paste(ifelse(is.na(row), "NA", row), collapse = ":"))
  length(unique(keys))
}

truth_for_dgp <- function(spec, seed, truth_n) {
  set.seed(seed)
  X <- simulate_ratings(truth_n, spec)
  rows <- list()
  weights <- unique(method_defs[, c("weight", "weight_label")])
  for (i in seq_len(nrow(weights))) {
    w <- weights$weight[i]
    wl <- weights$weight_label[i]
    fit <- misskappa::kappa(
      X,
      estimator = "ipw",
      weight = w,
      values = seq(0, spec$C - 1L)
    )
    est <- stats::coef(fit)
    rows[[length(rows) + 1L]] <- data.frame(
      dgp = spec$id,
      dgp_label = spec$label,
      dgp_family = spec$family,
      C = spec$C,
      R = spec$R,
      weight = w,
      weight_label = wl,
      coefficient = coef_names,
      truth = unname(est[coef_names]),
      truth_n = truth_n,
      truth_method = "complete-data Monte Carlo using kappa(estimator='ipw')",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

split_keys <- function(data, keys) {
  interaction(data[, keys], drop = TRUE, lex.order = TRUE)
}

summarize_replicates <- function(df) {
  keys <- c("dgp", "dgp_label", "dgp_family", "C", "R", "mechanism",
            "mechanism_family", "n", "method", "estimator", "weight_label",
            "coefficient")
  interval_names <- c("wald_z", "wald_t", "fisher_t", "asin_t")
  pieces <- lapply(split(df, split_keys(df, keys)), function(g) {
    ok <- is.finite(g$estimate)
    se_ok <- ok & is.finite(g$se) & g$se > 0
    err <- g$estimate - g$truth
    base <- data.frame(
      dgp = g$dgp[[1L]],
      dgp_label = g$dgp_label[[1L]],
      dgp_family = g$dgp_family[[1L]],
      C = g$C[[1L]],
      R = g$R[[1L]],
      mechanism = g$mechanism[[1L]],
      mechanism_family = g$mechanism_family[[1L]],
      n = g$n[[1L]],
      method = g$method[[1L]],
      estimator = g$estimator[[1L]],
      weight_label = g$weight_label[[1L]],
      coefficient = g$coefficient[[1L]],
      reps = length(unique(g$rep)),
      n_valid = sum(ok),
      failures = sum(!ok),
      truth = g$truth[[1L]],
      mean_estimate = if (any(ok)) mean(g$estimate[ok]) else NA_real_,
      bias = if (any(ok)) mean(err[ok]) else NA_real_,
      sd_estimate = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) else NA_real_,
      mc_se_bias = if (sum(ok) > 1L) stats::sd(g$estimate[ok]) / sqrt(sum(ok)) else NA_real_,
      mse = if (any(ok)) mean(err[ok]^2) else NA_real_,
      rmse = if (any(ok)) sqrt(mean(err[ok]^2)) else NA_real_,
      mean_se = if (any(se_ok)) mean(g$se[se_ok]) else NA_real_,
      se_over_sd = if (sum(se_ok) > 1L && stats::sd(g$estimate[se_ok]) > 1e-12) {
        mean(g$se[se_ok]) / stats::sd(g$estimate[se_ok])
      } else NA_real_,
      mean_elapsed_ms = mean(g$elapsed_ms, na.rm = TRUE),
      median_elapsed_ms = stats::median(g$elapsed_ms, na.rm = TRUE),
      mean_observed_fraction = mean(g$observed_fraction, na.rm = TRUE),
      mean_subjects_used = mean(g$subjects_used, na.rm = TRUE),
      mean_empty_rows = mean(g$empty_rows, na.rm = TRUE),
      mean_complete_rows = mean(g$complete_rows, na.rm = TRUE),
      mean_n_eff = mean(g$n_eff, na.rm = TRUE),
      min_pair_count_min = min(g$min_pair_count, na.rm = TRUE),
      mean_observed_patterns = mean(g$observed_patterns, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    do.call(rbind, lapply(interval_names, function(interval) {
      lo <- g[[paste0(interval, "_lower")]]
      hi <- g[[paste0(interval, "_upper")]]
      ci_ok <- se_ok & is.finite(lo) & is.finite(hi)
      covered <- ci_ok & lo <= g$truth & g$truth <= hi
      below <- ci_ok & hi < g$truth
      above <- ci_ok & lo > g$truth
      cbind(
        base,
        data.frame(
          interval = interval,
          interval_n = sum(ci_ok),
          coverage95 = if (any(ci_ok)) mean(covered[ci_ok]) else NA_real_,
          miss_below = if (any(ci_ok)) mean(below[ci_ok]) else NA_real_,
          miss_above = if (any(ci_ok)) mean(above[ci_ok]) else NA_real_,
          mean_ci_length = if (any(ci_ok)) mean(hi[ci_ok] - lo[ci_ok]) else NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }))
  })
  ans <- do.call(rbind, pieces)
  ans <- ans[order(ans$dgp, ans$mechanism, ans$n, ans$weight_label,
                   ans$estimator, ans$coefficient, ans$interval), ]
  rownames(ans) <- NULL
  ans
}

log_progress <- function(...) {
  if (opts$progress) message(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
}

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE)
  if (file.exists(path)) unlink(path)
  ok <- file.rename(tmp, path)
  if (!ok) stop("Failed to move temporary file into place: ", path, call. = FALSE)
}

checkpoint_file <- function(cell, grid, dir) {
  file.path(
    dir,
    sprintf(
      "cell-%03d-%s-%s-n%d-shard%02dof%02d.csv",
      cell,
      grid$dgp[[cell]],
      grid$mechanism[[cell]],
      grid$n[[cell]],
      opts$shard_index,
      opts$shard_count
    )
  )
}

method_is_applicable <- function(spec, method_row) {
  if (method_row$base_estimator %in% c("pairwise", "nt_fiml") &&
      method_row$weight != "quadratic") {
    return(FALSE)
  }
  TRUE
}

rep_ids_for_shard <- function(reps, shard_index, shard_count) {
  ids <- seq_len(reps)
  ids[((ids - 1L) %% shard_count) + 1L == shard_index]
}

mechanism_families <- c(
  complete = "sanity",
  mcar15 = "mcar",
  mcar30 = "mcar",
  mcar_rater30 = "mcar",
  mcar_subject30 = "mcar",
  mar_anchor15 = "mar",
  mar_anchor30 = "mar",
  mar_anchor_nonlinear30 = "mar",
  mar_shifted30 = "mar",
  mar_sequential30 = "mar",
  designed_bib3 = "planned",
  designed_forms4 = "planned",
  designed_anchor3 = "planned",
  designed_two_phase = "planned"
)

mechanism_descriptions <- c(
  complete = "No deletion; complete-data sanity row.",
  mcar15 = "Independent per-cell MCAR with 15 percent deletion.",
  mcar30 = "Independent per-cell MCAR with 30 percent deletion.",
  mcar_rater30 = "Rater-specific MCAR deletion rates with mean near 30 percent.",
  mcar_subject30 = "Subject-clustered MCAR deletion propensities with mean near 30 percent.",
  mar_anchor15 = "Clean anchor MAR: rater 1 always observed; later deletion depends linearly on rater 1 category, target 15 percent.",
  mar_anchor30 = "Clean anchor MAR: rater 1 always observed; later deletion depends linearly on rater 1 category, target 30 percent.",
  mar_anchor_nonlinear30 = "Clean nonlinear anchor MAR: later deletion is highest at extreme rater 1 categories, target 30 percent.",
  mar_shifted30 = "Clean shifted anchor MAR: anchor dependence plus rater-specific deletion shifts, target 30 percent on average.",
  mar_sequential30 = "Sequential MAR: later deletion depends on the previous observed rating or its missingness status, target 30 percent.",
  designed_bib3 = "Planned balanced incomplete block style: each subject receives exactly three raters when possible.",
  designed_forms4 = "Planned overlapping forms: fixed four-rater subsets cycle across subjects.",
  designed_anchor3 = "Planned anchor design: rater 1 plus two rotating additional raters.",
  designed_two_phase = "Planned two-phase design: anchor-three sparse phase plus every fifth subject fully rated."
)

truth_rows <- do.call(rbind, lapply(seq_along(opts$dgps), function(i) {
  spec <- dgp_defs[[opts$dgps[[i]]]]
  truth_for_dgp(spec, opts$seed_base + 900000L + i, opts$truth_n)
}))

grid <- expand.grid(
  dgp = opts$dgps,
  mechanism = opts$mechanisms,
  n = opts$n_grid,
  stringsAsFactors = FALSE
)
grid <- grid[order(grid$dgp, grid$mechanism, grid$n), , drop = FALSE]
grid$cell <- seq_len(nrow(grid))

mechanisms_out <- data.frame(
  mechanism = valid_mechanisms,
  family = unname(mechanism_families[valid_mechanisms]),
  description = unname(mechanism_descriptions[valid_mechanisms]),
  stringsAsFactors = FALSE
)

method_rows <- method_defs

detail_string <- function(spec) {
  keys <- names(spec)[names(spec) %in% c("target_quadratic", "knowledge_prob", "latent_corr")]
  if (!length(keys)) return("")
  paste(keys, signif(as.numeric(unlist(spec[keys])), 4), collapse = ";")
}

dgp_rows <- do.call(rbind, lapply(dgp_defs[opts$dgps], function(spec) {
  data.frame(
    dgp = spec$id,
    dgp_label = spec$label,
    family = spec$family,
    C = spec$C,
    R = spec$R,
    support_size = spec$C^spec$R,
    category_prob = if (!is.null(spec$p)) paste(signif(spec$p, 4), collapse = ",") else "",
    details = detail_string(spec),
    stringsAsFactors = FALSE
  )
}))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(opts$out_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_atomic(truth_rows, file.path(opts$out_dir, "truth.csv"))
write_csv_atomic(grid[, c("cell", "dgp", "mechanism", "n")],
                 file.path(opts$out_dir, "cell_plan.csv"))
write_csv_atomic(mechanisms_out, file.path(opts$out_dir, "mechanisms.csv"))
write_csv_atomic(method_rows, file.path(opts$out_dir, "methods.csv"))
write_csv_atomic(dgp_rows, file.path(opts$out_dir, "dgps.csv"))

checkpoint_files <- character(nrow(grid))
t0 <- Sys.time()
shard_rep_ids <- rep_ids_for_shard(opts$reps, opts$shard_index, opts$shard_count)
for (cell in seq_len(nrow(grid))) {
  spec <- dgp_defs[[grid$dgp[[cell]]]]
  checkpoint_files[[cell]] <- checkpoint_file(cell, grid, checkpoint_dir)
  if (opts$resume && file.exists(checkpoint_files[[cell]])) {
    log_progress("cell %d/%d: reusing %s", cell, nrow(grid), basename(checkpoint_files[[cell]]))
    next
  }

  cell_rows <- list()
  log_progress(
    "cell %d/%d: %s %s n=%d reps=%d shard=%d/%d",
    cell, nrow(grid), grid$dgp[[cell]], grid$mechanism[[cell]],
    grid$n[[cell]], length(shard_rep_ids), opts$shard_index, opts$shard_count
  )
  for (rep_id in shard_rep_ids) {
    seed <- opts$seed_base + 1000000L * cell + rep_id
    set.seed(seed)
    X_star <- simulate_ratings(grid$n[[cell]], spec)
    X <- apply_missing(X_star, spec, grid$mechanism[[cell]])
    observed_fraction <- mean(!is.na(X))
    empty_rows <- sum(rowSums(!is.na(X)) == 0L)
    keep <- rowSums(!is.na(X)) > 0L
    X <- X[keep, , drop = FALSE]
    subjects_used <- nrow(X)
    complete_rows <- sum(stats::complete.cases(X))
    min_pair <- min_pair_count(X)
    observed_patterns <- observed_pattern_count(X)

    for (m in seq_len(nrow(method_defs))) {
      method_row <- method_defs[m, ]
      if (!method_is_applicable(spec, method_row)) next
      fit <- fit_one(X, spec, method_row)
      truth_key <- truth_rows$dgp == spec$id &
        truth_rows$weight_label == method_row$weight_label
      truth_now <- truth_rows[truth_key, ]
      fit$truth <- truth_now$truth[match(fit$coefficient, truth_now$coefficient)]

      cell_rows[[length(cell_rows) + 1L]] <- data.frame(
        dgp = spec$id,
        dgp_label = spec$label,
        dgp_family = spec$family,
        C = spec$C,
        R = spec$R,
        mechanism = grid$mechanism[[cell]],
        mechanism_family = unname(mechanism_families[[grid$mechanism[[cell]]]]),
        n = grid$n[[cell]],
        rep = rep_id,
        seed = seed,
        method = method_row$method,
        estimator = method_row$estimator,
        weight = method_row$weight,
        weight_label = method_row$weight_label,
        fit,
        observed_fraction = observed_fraction,
        subjects_used = subjects_used,
        empty_rows = empty_rows,
        complete_rows = complete_rows,
        min_pair_count = min_pair,
        observed_patterns = observed_patterns,
        stringsAsFactors = FALSE
      )
    }
  }
  cell_replicates <- if (length(cell_rows)) do.call(rbind, cell_rows) else data.frame()
  write_csv_atomic(cell_replicates, checkpoint_files[[cell]])
  log_progress("cell %d/%d: wrote %s (%d rows)",
               cell, nrow(grid), basename(checkpoint_files[[cell]]), nrow(cell_replicates))
}

missing_checkpoints <- checkpoint_files[!file.exists(checkpoint_files)]
if (length(missing_checkpoints)) {
  stop("Missing checkpoint(s): ", paste(missing_checkpoints, collapse = ", "), call. = FALSE)
}
replicates <- do.call(rbind, lapply(checkpoint_files, function(path) {
  read.csv(path, stringsAsFactors = FALSE)
}))
summary <- summarize_replicates(replicates)

write_csv_atomic(replicates, file.path(opts$out_dir, "replicates.csv"))
write_csv_atomic(summary, file.path(opts$out_dir, "summary.csv"))

cmd <- paste(commandArgs(FALSE), collapse = " ")
metadata <- data.frame(
  key = c("generated_at", "script", "command", "mode", "reps", "reps_in_shard",
          "n_grid", "truth_n", "dgps", "mechanisms", "methods", "seed_base",
          "shard_index", "shard_count", "resume", "checkpoint_dir",
          "elapsed_seconds", "misskappa_version", "r_version"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    normalizePath(file.path(script_dir, "run_experiment.R"), mustWork = FALSE),
    cmd,
    opts$mode,
    as.character(opts$reps),
    as.character(length(shard_rep_ids)),
    paste(opts$n_grid, collapse = ","),
    as.character(opts$truth_n),
    paste(opts$dgps, collapse = ","),
    paste(opts$mechanisms, collapse = ","),
    paste(opts$methods, collapse = ","),
    as.character(opts$seed_base),
    as.character(opts$shard_index),
    as.character(opts$shard_count),
    as.character(opts$resume),
    normalizePath(checkpoint_dir, mustWork = FALSE),
    as.character(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)),
    as.character(utils::packageVersion("misskappa")),
    R.version.string
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(metadata, file.path(opts$out_dir, "metadata.csv"))

cat("Wrote:\n")
for (name in c("truth.csv", "cell_plan.csv", "replicates.csv", "summary.csv",
               "mechanisms.csv", "methods.csv", "dgps.csv", "metadata.csv")) {
  cat("  ", normalizePath(file.path(opts$out_dir, name), mustWork = FALSE), "\n", sep = "")
}
cat("  ", normalizePath(checkpoint_dir, mustWork = FALSE), "/\n", sep = "")

headline <- summary[summary$coefficient == "Conger" &
                      summary$interval == "wald_t" &
                      summary$weight_label == "quadratic",
                    c("dgp", "mechanism", "method", "bias", "rmse",
                      "coverage95", "n_valid", "failures", "mean_elapsed_ms")]
headline <- headline[order(headline$dgp, headline$mechanism, headline$method), ]
if (nrow(headline)) {
  cat("\nConger quadratic smoke rows, Wald-t interval:\n")
  print(utils::head(headline, 24L), row.names = FALSE, digits = 3)
}
