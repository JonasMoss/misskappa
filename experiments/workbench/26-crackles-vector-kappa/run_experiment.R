#!/usr/bin/env Rscript
#
# 26-crackles-vector-kappa
#
# CRACKLES pilot and synthetic verification for the internal
# component-separable vector kappa estimator and the full-weight quadratic
# covariance route.

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}
if (has_flag("--help") || has_flag("-h")) {
  cat("Usage: Rscript run_experiment.R [options]\n",
      " --mask-prop P   Synthetic component-missing fraction (default 0.25).\n",
      " --reps N        Number of masking replicates (default 25; smoke 3).\n",
      " --quad-n N      Normal quadratic simulation n (default 250; smoke 160).\n",
      " --quad-reps N   Normal quadratic simulation reps (default --reps).\n",
      " --seed-base N   Deterministic seed base (default 26000).\n",
      " --smoke         Cheap run with 3 masking replicates and smaller synthetic n.\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

mask_prop <- get_val("--mask-prop", 0.25, as.numeric)
reps <- get_val("--reps", if (has_flag("--smoke")) 3L else 25L, as.integer)
quad_reps <- get_val("--quad-reps", reps, as.integer)
seed_base <- get_val("--seed-base", 26000L, as.integer)
synthetic_n <- if (has_flag("--smoke")) 80L else 500L
quadratic_n <- get_val("--quad-n", if (has_flag("--smoke")) 160L else 250L,
                       as.integer)

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1L])
script_dir <- if (length(script_file) && !is.na(script_file)) {
  dirname(normalizePath(script_file))
} else {
  getwd()
}
results_dir <- file.path(script_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

find_repo_root <- function(start) {
  cur <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "AGENTS.md")) &&
        dir.exists(file.path(cur, "r-package"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not find repository root from ", start, call. = FALSE)
    }
    cur <- parent
  }
}

repo_root <- find_repo_root(script_dir)
data_path <- file.path(repo_root, "dev", "dat", "CRACKLES.rda")
if (!file.exists(data_path)) {
  stop("Missing CRACKLES data at ", data_path, call. = FALSE)
}

pkg_dir <- file.path(repo_root, "r-package")
if (requireNamespace("pkgload", quietly = TRUE) &&
    file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  pkgload::load_all(pkg_dir, export_all = FALSE, quiet = TRUE)
} else {
  library(misskappa)
}
vector_kappa <- getFromNamespace("kappa_vector", "misskappa")
vector_quadratic <- getFromNamespace("kappa_vector_quadratic", "misskappa")
kvq_grad <- getFromNamespace(".kvq_grad", "misskappa")

load(data_path)
d <- CRACKLES
rating_cols <- setdiff(names(d), c("patient", "UP", "LO"))
groups <- unique(sub("[0-9]+$", "", rating_cols))
observer_classes <- data.frame(
  observer_code = c("ALL", "EXP", "NOR", "RUS", "WAL", "NLD", "PUL", "STU"),
  observer_class = c(
    "All observer classes",
    "International lung-sound experts / researchers",
    "General practitioners from Norway",
    "General practitioners from Russia",
    "General practitioners from Wales",
    "General practitioners from The Netherlands",
    "Pulmonologists, University Hospital of North Norway",
    "Sixth-year medical students in Tromso"
  ),
  vanbelle_label = c("All", "EXP", "NOR", "RUS", "WAL", "NLD", "PLN", "STU"),
  stringsAsFactors = FALSE
)

make_site_layout <- function(d) {
  split_idx <- split(seq_len(nrow(d)), d$patient)
  n_by_patient <- vapply(split_idx, length, integer(1))
  if (length(unique(n_by_patient)) != 1L) {
    stop("Patients do not all have the same number of rows.", call. = FALSE)
  }
  site_ord <- ave(seq_len(nrow(d)), d$patient, FUN = seq_along)
  location_code <- ifelse(d$UP == 1L, "U", ifelse(d$LO == 1L, "L", "A"))
  location <- ifelse(d$UP == 1L, "upper posterior",
                     ifelse(d$LO == 1L, "lower posterior", "anterior"))
  side_index <- ave(location_code, d$patient, location_code, FUN = seq_along)
  data.frame(
    row = seq_len(nrow(d)),
    patient = d$patient,
    site = site_ord,
    location_code = location_code,
    location = location,
    side_index = side_index,
    site_label = paste0(location_code, side_index),
    stringsAsFactors = FALSE
  )
}

pivot_array <- function(d, cols, site_labels) {
  patients <- sort(unique(d$patient))
  split_idx <- split(seq_len(nrow(d)), d$patient)
  features <- length(split_idx[[1L]])
  R <- length(cols)
  X <- array(NA_real_, dim = c(length(patients), R, features),
             dimnames = list(patient = patients, rater = cols,
                             site = site_labels))
  for (i in seq_along(patients)) {
    rows <- split_idx[[as.character(patients[i])]]
    X[i, , ] <- t(as.matrix(d[rows, cols]))
  }
  storage.mode(X) <- "double"
  X
}

fit_vector <- function(X, method, loss = "hamming", feature_weights = NULL) {
  fit <- vector_kappa(X, method = method, loss = loss,
                      feature_weights = feature_weights)
  V <- vcov(fit)
  data.frame(
    method = method,
    loss = loss,
    kappa_C = unname(fit$estimates["Conger"]),
    kappa_F = unname(fit$estimates["Fleiss"]),
    se_C = sqrt(V["Conger", "Conger"]),
    se_F = sqrt(V["Fleiss", "Fleiss"]),
    stringsAsFactors = FALSE
  )
}

fit_grid <- function(X, losses = "hamming", methods = c("pairwise", "ipw"),
                     feature_weights = NULL) {
  rows <- list()
  for (loss in losses) {
    for (method in methods) {
      rows[[length(rows) + 1L]] <- fit_vector(X, method, loss, feature_weights)
    }
  }
  do.call(rbind, rows)
}

make_site_W <- function(features) {
  W <- diag(features)
  if (features == 6L) {
    for (idx in list(c(1L, 2L), c(3L, 4L), c(5L, 6L))) {
      W[idx[1L], idx[2L]] <- W[idx[2L], idx[1L]] <- 0.25
    }
  }
  W
}

fit_quadratic <- function(X, method, W, em_options = list(),
                          loss = "quadratic_full_W") {
  fit <- vector_quadratic(X, method = method, W = W, em_options = em_options)
  V <- vcov(fit)
  data.frame(
    method = method,
    loss = loss,
    kappa_C = unname(fit$estimates["Conger"]),
    kappa_F = unname(fit$estimates["Fleiss"]),
    se_C = sqrt(V["Conger", "Conger"]),
    se_F = sqrt(V["Fleiss", "Fleiss"]),
    n_used = dim(X)[1L],
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

safe_fit_quadratic <- function(X, method, W, em_options = list(),
                               loss = "quadratic_full_W") {
  tryCatch(
    fit_quadratic(X, method, W, em_options = em_options, loss = loss),
    error = function(e) {
      data.frame(
        method = method,
        loss = loss,
        kappa_C = NA_real_,
        kappa_F = NA_real_,
        se_C = NA_real_,
        se_F = NA_real_,
        n_used = dim(X)[1L],
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

fit_quadratic_methods <- function(X, W, methods = "pairwise",
                                  em_options = list()) {
  rows <- lapply(methods, function(method) {
    safe_fit_quadratic(X, method, W, em_options = em_options)
  })
  do.call(rbind, rows)
}

mask_mcar <- function(X, prop) {
  Xm <- X
  drop <- stats::runif(length(Xm)) < prop
  Xm[drop] <- NA_real_
  Xm
}

fit_masked_grid <- function(X, prop, losses = "hamming",
                            feature_weights = NULL, max_tries = 100L) {
  for (try in seq_len(max_tries)) {
    Xm <- mask_mcar(X, prop)
    fits <- try(fit_grid(Xm, losses = losses, feature_weights = feature_weights),
                silent = TRUE)
    if (!inherits(fits, "try-error")) {
      fits$missing_fraction <- mean(is.na(Xm))
      return(fits)
    }
  }
  stop("Could not generate a non-singular masked data set.", call. = FALSE)
}

summarise_group <- function(g) {
  vals <- c("kappa_C", "kappa_F", "err_C", "err_F", "se_C", "se_F",
            "missing_fraction")
  out <- g[1L, intersect(c("analysis", "analysis_label", "observer_code",
                           "observer_class", "vanbelle_label", "method",
                           "loss", "R", "features", "dgp"), names(g))]
  for (v in vals) {
    if (!v %in% names(g)) next
    x <- as.numeric(g[[v]])
    ok <- !is.na(x)
    out[[paste0(v, "_n_valid")]] <- sum(ok)
    out[[paste0(v, "_mean")]] <- if (any(ok)) mean(x[ok]) else NA_real_
    out[[paste0(v, "_sd")]] <- if (sum(ok) > 1L) stats::sd(x[ok]) else NA_real_
    out[[paste0(v, "_q05")]] <- if (any(ok)) stats::quantile(x[ok], 0.05) else NA_real_
    out[[paste0(v, "_q95")]] <- if (any(ok)) stats::quantile(x[ok], 0.95) else NA_real_
  }
  out
}

summarise_quadratic_sim <- function(g, truth) {
  coef_info <- data.frame(
    suffix = c("C", "F"),
    coefficient = c("Conger", "Fleiss"),
    truth = c(unname(truth["Conger"]), unname(truth["Fleiss"])),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(coef_info)), function(i) {
    s <- coef_info$suffix[i]
    est <- as.numeric(g[[paste0("kappa_", s)]])
    se <- as.numeric(g[[paste0("se_", s)]])
    ok <- is.finite(est)
    ok_se <- ok & is.finite(se)
    err <- est - coef_info$truth[i]
    cover <- abs(err) <= 1.96 * se
    data.frame(
      mechanism = g$mechanism[1L],
      method = g$method[1L],
      loss = g$loss[1L],
      coefficient = coef_info$coefficient[i],
      truth = coef_info$truth[i],
      valid_reps = sum(ok),
      bias = if (any(ok)) mean(err[ok]) else NA_real_,
      rmse = if (any(ok)) sqrt(mean(err[ok]^2)) else NA_real_,
      sd_est = if (sum(ok) > 1L) stats::sd(est[ok]) else NA_real_,
      mean_se = if (any(ok_se)) mean(se[ok_se]) else NA_real_,
      coverage95 = if (any(ok_se)) mean(cover[ok_se]) else NA_real_,
      missing_fraction_mean = mean(g$missing_fraction),
      n_used_mean = mean(g$n_used[ok], na.rm = TRUE),
      error_reps = sum(!is.na(g$error)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

site_layout <- make_site_layout(d)
site_summary <- unique(site_layout[, c("site", "location_code", "location",
                                       "side_index", "site_label")])
write.csv(site_summary, file.path(results_dir, "site_layout.csv"), row.names = FALSE)

analysis_sets <- c(list(pooled = rating_cols),
                   setNames(lapply(groups, function(g) rating_cols[startsWith(rating_cols, g)]),
                            paste0("panel_", groups)))
analysis_info <- data.frame(
  analysis = names(analysis_sets),
  observer_code = c("ALL", groups),
  stringsAsFactors = FALSE
)
analysis_info <- merge(analysis_info, observer_classes, by = "observer_code",
                       sort = FALSE)
analysis_info$analysis_label <- ifelse(
  analysis_info$observer_code == "ALL",
  "All observer classes (descriptive)",
  paste0(analysis_info$observer_code, " - ", analysis_info$observer_class)
)
analysis_info <- analysis_info[match(names(analysis_sets), analysis_info$analysis), ]

complete_rows <- list()
masked_rows <- list()
quadratic_complete_rows <- list()
set.seed(seed_base)
for (nm in names(analysis_sets)) {
  cols <- analysis_sets[[nm]]
  info <- analysis_info[analysis_info$analysis == nm, ]
  X <- pivot_array(d, cols, site_summary$site_label)
  R <- dim(X)[2L]
  features <- dim(X)[3L]
  feature_weights <- rep(1, features)

  full <- fit_grid(X, losses = "hamming", feature_weights = feature_weights)
  complete_rows[[length(complete_rows) + 1L]] <- cbind(
    analysis = nm, analysis_label = info$analysis_label,
    observer_code = info$observer_code, observer_class = info$observer_class,
    vanbelle_label = info$vanbelle_label,
    n_patients = dim(X)[1L], R = R, features = features,
    missing_fraction = mean(is.na(X)), full)

  W_site <- make_site_W(features)
  qfull <- fit_quadratic_methods(X, W = W_site, methods = "pairwise")
  quadratic_complete_rows[[length(quadratic_complete_rows) + 1L]] <- cbind(
    analysis = nm, analysis_label = info$analysis_label,
    observer_code = info$observer_code, observer_class = info$observer_class,
    vanbelle_label = info$vanbelle_label,
    n_patients = dim(X)[1L], R = R, features = features,
    missing_fraction = mean(is.na(X)),
    W_rank = qr(W_site)$rank, W_trace = sum(diag(W_site)), qfull)

  truth <- full[full$method == "pairwise" & full$loss == "hamming",
                c("kappa_C", "kappa_F")]
  names(truth) <- c("truth_C", "truth_F")
  for (rep in seq_len(reps)) {
    fits <- fit_masked_grid(X, mask_prop, losses = "hamming",
                            feature_weights = feature_weights)
    fits$err_C <- fits$kappa_C - truth$truth_C
    fits$err_F <- fits$kappa_F - truth$truth_F
    masked_rows[[length(masked_rows) + 1L]] <- cbind(
      analysis = nm, analysis_label = info$analysis_label,
      observer_code = info$observer_code, observer_class = info$observer_class,
      vanbelle_label = info$vanbelle_label,
      rep = rep, n_patients = dim(X)[1L], R = R, features = features,
      fits)
  }
}

simulate_binary_vector <- function(n, R, features, seed) {
  set.seed(seed)
  feature_prob <- seq(0.25, 0.75, length.out = features)
  truth <- matrix(stats::rbinom(n * features, 1, rep(feature_prob, each = n)),
                  nrow = n, ncol = features)
  skill <- seq(0.92, 0.68, length.out = R)
  X <- array(NA_real_, dim = c(n, R, features))
  for (j in seq_len(R)) {
    flip <- matrix(stats::rbinom(n * features, 1, 1 - skill[j]),
                   nrow = n, ncol = features)
    X[, j, ] <- abs(truth - flip)
  }
  X
}

mask_heterogeneous <- function(X, miss_prob) {
  Xm <- X
  for (j in seq_len(dim(X)[2L])) {
    for (l in seq_len(dim(X)[3L])) {
      drop <- stats::runif(dim(X)[1L]) < miss_prob[j, l]
      Xm[drop, j, l] <- NA_real_
    }
  }
  Xm
}

synthetic_complete <- simulate_binary_vector(synthetic_n, 4L, 3L, seed_base + 7000L)
synthetic_weights <- c(1, 2, 4)
synthetic_truth <- fit_grid(synthetic_complete, losses = c("hamming", "rms"),
                            feature_weights = synthetic_weights)
synthetic_truth$dgp <- "heterogeneous_mcar"
synthetic_truth$rep <- 0L
synthetic_truth$missing_fraction <- 0
synthetic_truth$err_C <- 0
synthetic_truth$err_F <- 0

miss_prob <- matrix(c(
  0.10, 0.45, 0.70,
  0.15, 0.35, 0.65,
  0.20, 0.30, 0.55,
  0.25, 0.25, 0.45
), nrow = 4L, byrow = TRUE)
synthetic_rows <- list()
set.seed(seed_base + 8000L)
for (rep in seq_len(reps)) {
  Xm <- mask_heterogeneous(synthetic_complete, miss_prob)
  fits <- fit_grid(Xm, losses = c("hamming", "rms"),
                   feature_weights = synthetic_weights)
  for (loss in unique(fits$loss)) {
    truth <- synthetic_truth[synthetic_truth$method == "pairwise" &
                               synthetic_truth$loss == loss,
                             c("kappa_C", "kappa_F")]
    idx <- fits$loss == loss
    fits$err_C[idx] <- fits$kappa_C[idx] - truth$kappa_C
    fits$err_F[idx] <- fits$kappa_F[idx] - truth$kappa_F
  }
  synthetic_rows[[length(synthetic_rows) + 1L]] <- cbind(
    dgp = "heterogeneous_mcar", rep = rep, n_patients = dim(Xm)[1L],
    R = dim(Xm)[2L], features = dim(Xm)[3L],
    missing_fraction = mean(is.na(Xm)), fits)
}

make_quadratic_dgp <- function(R = 4L, features = 3L) {
  Lw <- matrix(c(
    1.00, 0.00, 0.00,
    0.25, 1.05, 0.00,
   -0.10, 0.35, 1.15
  ), nrow = features, byrow = TRUE)
  W <- crossprod(Lw)

  Lf <- matrix(c(
    1.00, 0.00, 0.00,
    0.35, 0.90, 0.00,
    0.15, 0.25, 0.80
  ), nrow = features, byrow = TRUE)
  feature_cov <- crossprod(Lf)
  rater_idx <- seq_len(R)
  rater_cov <- 0.58 ^ abs(outer(rater_idx, rater_idx, "-"))
  residual <- diag(c(0.30, 0.22, 0.36), features)
  Sigma <- kronecker(rater_cov, feature_cov) + kronecker(diag(R), residual)

  feature_level <- c(-0.45, 0.05, 0.55)
  rater_shift <- c(-0.25, -0.05, 0.10, 0.24)
  mu_mat <- outer(rater_shift, rep(1, features)) +
    outer(rep(1, R), feature_level)
  mu <- as.vector(t(mu_mat))
  list(R = R, features = features, W = W, mu = mu, Sigma = Sigma)
}

simulate_quadratic_array <- function(n, dgp, seed) {
  set.seed(seed)
  q <- length(dgp$mu)
  Xflat <- matrix(stats::rnorm(n * q), nrow = n) %*% chol(dgp$Sigma)
  Xflat <- sweep(Xflat, 2L, dgp$mu, "+")
  X <- array(NA_real_, dim = c(n, dgp$R, dgp$features))
  for (j in seq_len(dgp$R)) {
    for (l in seq_len(dgp$features)) {
      X[, j, l] <- Xflat[, (j - 1L) * dgp$features + l]
    }
  }
  X
}

std_score <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

mask_quadratic_mcar <- function(X) {
  Xm <- X
  n <- dim(X)[1L]; R <- dim(X)[2L]; features <- dim(X)[3L]
  component_prob <- matrix(c(
    0.02, 0.06, 0.12,
    0.04, 0.08, 0.14,
    0.06, 0.10, 0.16,
    0.08, 0.12, 0.20
  ), nrow = R, byrow = TRUE)
  profile_prob <- c(0.00, 0.02, 0.03, 0.05)
  profile_drop <- matrix(FALSE, n, R)
  for (j in seq_len(R)) {
    profile_drop[, j] <- stats::runif(n) < profile_prob[j]
  }
  for (j in seq_len(R)) {
    for (l in seq_len(features)) {
      drop <- profile_drop[, j] |
        stats::runif(n) < component_prob[j, l]
      Xm[drop, j, l] <- NA_real_
    }
  }
  Xm
}

mask_quadratic_mar <- function(X) {
  Xm <- X
  n <- dim(X)[1L]; R <- dim(X)[2L]; features <- dim(X)[3L]
  anchor_global <- std_score(rowMeans(X[, 1L, ]))
  anchor_feature <- apply(X[, 1L, ], 2L, std_score)
  component_base <- matrix(c(
    0.00, 0.00, 0.00,
    0.05, 0.10, 0.18,
    0.08, 0.14, 0.23,
    0.10, 0.18, 0.30
  ), nrow = R, byrow = TRUE)
  profile_base <- c(0.00, 0.02, 0.04, 0.06)
  profile_drop <- matrix(FALSE, n, R)
  for (j in 2L:R) {
    eta_profile <- stats::qlogis(profile_base[j]) +
      0.70 * anchor_global + 0.12 * (j - 2L)
    profile_drop[, j] <- stats::runif(n) < stats::plogis(eta_profile)
  }
  for (j in 2L:R) {
    for (l in seq_len(features)) {
      eta_component <- stats::qlogis(component_base[j, l]) +
        0.65 * anchor_global + 0.40 * anchor_feature[, l] +
        0.10 * (j - 2L) - 0.08 * (l - 2L)
      drop <- profile_drop[, j] |
        stats::runif(n) < stats::plogis(eta_component)
      Xm[drop, j, l] <- NA_real_
    }
  }
  Xm
}

fit_quadratic_sim_method <- function(X, method, W) {
  Xfit <- X
  backend <- method
  if (method == "listwise") {
    keep <- apply(is.finite(X), 1L, all)
    Xfit <- X[keep, , , drop = FALSE]
    backend <- "pairwise"
  }
  if (dim(Xfit)[1L] < 10L) {
    stop("fewer than 10 usable subjects.", call. = FALSE)
  }
  fit <- safe_fit_quadratic(
    Xfit, backend, W,
    em_options = list(tol = 1e-7, max_iter = 2000L, fd_h = 1e-4)
  )
  fit$method <- method
  fit$n_used <- dim(Xfit)[1L]
  fit
}

quadratic_dgp <- make_quadratic_dgp()
quadratic_truth <- kvq_grad(
  quadratic_dgp$mu, quadratic_dgp$Sigma,
  quadratic_dgp$R, quadratic_dgp$features, quadratic_dgp$W
)$estimates
quadratic_truth_df <- data.frame(
  coefficient = names(quadratic_truth),
  truth = unname(quadratic_truth),
  R = quadratic_dgp$R,
  features = quadratic_dgp$features,
  stringsAsFactors = FALSE
)

quadratic_sim_rows <- list()
mechanisms <- list(mcar_layered = mask_quadratic_mcar,
                   mar_anchor = mask_quadratic_mar)
for (rep in seq_len(quad_reps)) {
  Xq <- simulate_quadratic_array(
    quadratic_n, quadratic_dgp, seed_base + 9000L + rep
  )
  for (mech in names(mechanisms)) {
    Xmiss <- mechanisms[[mech]](Xq)
    for (method in c("listwise", "pairwise", "nt_fiml")) {
      fit <- tryCatch(
        fit_quadratic_sim_method(Xmiss, method, quadratic_dgp$W),
        error = function(e) {
          data.frame(
            method = method,
            loss = "quadratic_full_W",
            kappa_C = NA_real_,
            kappa_F = NA_real_,
            se_C = NA_real_,
            se_F = NA_real_,
            n_used = NA_integer_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
      fit$err_C <- fit$kappa_C - unname(quadratic_truth["Conger"])
      fit$err_F <- fit$kappa_F - unname(quadratic_truth["Fleiss"])
      quadratic_sim_rows[[length(quadratic_sim_rows) + 1L]] <- cbind(
        mechanism = mech, rep = rep, n_patients = dim(Xmiss)[1L],
        R = dim(Xmiss)[2L], features = dim(Xmiss)[3L],
        missing_fraction = mean(is.na(Xmiss)), fit)
    }
  }
}

complete_df <- do.call(rbind, complete_rows)
masked_df <- do.call(rbind, masked_rows)
quadratic_complete_df <- do.call(rbind, quadratic_complete_rows)
synthetic_df <- rbind(
  cbind(n_patients = dim(synthetic_complete)[1L], R = dim(synthetic_complete)[2L],
        features = dim(synthetic_complete)[3L], synthetic_truth),
  do.call(rbind, synthetic_rows)
)
quadratic_sim_df <- do.call(rbind, quadratic_sim_rows)

write.csv(complete_df, file.path(results_dir, "complete_estimates.csv"),
          row.names = FALSE)
write.csv(masked_df, file.path(results_dir, "masked_estimates.csv"),
          row.names = FALSE)
write.csv(synthetic_df, file.path(results_dir, "synthetic_estimates.csv"),
          row.names = FALSE)
write.csv(quadratic_complete_df,
          file.path(results_dir, "quadratic_complete_estimates.csv"),
          row.names = FALSE)
write.csv(quadratic_sim_df,
          file.path(results_dir, "quadratic_sim_estimates.csv"),
          row.names = FALSE)
write.csv(quadratic_truth_df[, c("coefficient", "truth", "R", "features")],
          file.path(results_dir, "quadratic_sim_truth.csv"),
          row.names = FALSE)

masked_keys <- unique(masked_df[, c("analysis", "method", "loss", "R", "features")])
masked_summary <- do.call(rbind, lapply(seq_len(nrow(masked_keys)), function(i) {
  keep <- masked_df$analysis == masked_keys$analysis[i] &
    masked_df$method == masked_keys$method[i] &
    masked_df$loss == masked_keys$loss[i]
  summarise_group(masked_df[keep, ])
}))
write.csv(masked_summary, file.path(results_dir, "masked_summary.csv"),
          row.names = FALSE)

synthetic_keys <- unique(synthetic_df[synthetic_df$rep > 0,
                                      c("dgp", "method", "loss", "R", "features")])
synthetic_summary <- do.call(rbind, lapply(seq_len(nrow(synthetic_keys)), function(i) {
  keep <- synthetic_df$rep > 0 &
    synthetic_df$dgp == synthetic_keys$dgp[i] &
    synthetic_df$method == synthetic_keys$method[i] &
    synthetic_df$loss == synthetic_keys$loss[i]
  summarise_group(synthetic_df[keep, ])
}))
write.csv(synthetic_summary, file.path(results_dir, "synthetic_summary.csv"),
          row.names = FALSE)

quadratic_sim_keys <- unique(quadratic_sim_df[
  , c("mechanism", "method", "loss", "R", "features")
])
quadratic_sim_summary <- do.call(rbind, lapply(seq_len(nrow(quadratic_sim_keys)),
                                               function(i) {
  keep <- quadratic_sim_df$mechanism == quadratic_sim_keys$mechanism[i] &
    quadratic_sim_df$method == quadratic_sim_keys$method[i] &
    quadratic_sim_df$loss == quadratic_sim_keys$loss[i]
  summarise_quadratic_sim(quadratic_sim_df[keep, ], quadratic_truth)
}))
write.csv(quadratic_sim_summary,
          file.path(results_dir, "quadratic_sim_summary.csv"),
          row.names = FALSE)

metadata <- data.frame(
  key = c("mask_prop", "reps", "seed_base", "synthetic_n", "n_rows",
          "quadratic_n", "quad_reps", "n_patients", "n_rating_columns", "R_version",
          "misskappa_version"),
  value = c(mask_prop, reps, seed_base, synthetic_n, nrow(d),
            quadratic_n, quad_reps, length(unique(d$patient)), length(rating_cols),
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(utils::packageVersion("misskappa"))),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(results_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n",
    " ", file.path(results_dir, "site_layout.csv"), "\n",
    " ", file.path(results_dir, "complete_estimates.csv"), "\n",
    " ", file.path(results_dir, "masked_estimates.csv"), "\n",
    " ", file.path(results_dir, "masked_summary.csv"), "\n",
    " ", file.path(results_dir, "synthetic_estimates.csv"), "\n",
    " ", file.path(results_dir, "synthetic_summary.csv"), "\n",
    " ", file.path(results_dir, "quadratic_complete_estimates.csv"), "\n",
    " ", file.path(results_dir, "quadratic_sim_estimates.csv"), "\n",
    " ", file.path(results_dir, "quadratic_sim_summary.csv"), "\n",
    " ", file.path(results_dir, "quadratic_sim_truth.csv"), "\n",
    " ", file.path(results_dir, "metadata.csv"), "\n", sep = "")
