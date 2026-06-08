#!/usr/bin/env Rscript
#
# 26-crackles-vector-kappa
#
# CRACKLES pilot and synthetic verification for the internal
# component-separable vector kappa estimator.

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
      " --seed-base N   Deterministic seed base (default 26000).\n",
      " --smoke         Cheap run with 3 masking replicates and smaller synthetic n.\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

mask_prop <- get_val("--mask-prop", 0.25, as.numeric)
reps <- get_val("--reps", if (has_flag("--smoke")) 3L else 25L, as.integer)
seed_base <- get_val("--seed-base", 26000L, as.integer)
synthetic_n <- if (has_flag("--smoke")) 80L else 500L

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1L])
script_dir <- if (length(script_file) && !is.na(script_file)) {
  dirname(normalizePath(script_file))
} else {
  getwd()
}
results_dir <- file.path(script_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

repo_root <- normalizePath(file.path(script_dir, "..", ".."))
data_path <- file.path(repo_root, "dev", "dat", "CRACKLES.rda")
if (!file.exists(data_path)) {
  stop("Missing CRACKLES data at ", data_path, call. = FALSE)
}

library(misskappa)
vector_kappa <- getFromNamespace("kappa_vector", "misskappa")

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

complete_df <- do.call(rbind, complete_rows)
masked_df <- do.call(rbind, masked_rows)
synthetic_df <- rbind(
  cbind(n_patients = dim(synthetic_complete)[1L], R = dim(synthetic_complete)[2L],
        features = dim(synthetic_complete)[3L], synthetic_truth),
  do.call(rbind, synthetic_rows)
)

write.csv(complete_df, file.path(results_dir, "complete_estimates.csv"),
          row.names = FALSE)
write.csv(masked_df, file.path(results_dir, "masked_estimates.csv"),
          row.names = FALSE)
write.csv(synthetic_df, file.path(results_dir, "synthetic_estimates.csv"),
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

metadata <- data.frame(
  key = c("mask_prop", "reps", "seed_base", "synthetic_n", "n_rows",
          "n_patients", "n_rating_columns", "R_version",
          "misskappa_version"),
  value = c(mask_prop, reps, seed_base, synthetic_n, nrow(d),
            length(unique(d$patient)), length(rating_cols),
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
    " ", file.path(results_dir, "metadata.csv"), "\n", sep = "")
