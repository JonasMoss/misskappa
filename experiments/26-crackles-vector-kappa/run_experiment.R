#!/usr/bin/env Rscript
#
# 26-crackles-vector-kappa
#
# Pure-R pilot for a vector-valued Fleiss/Conger agreement coefficient on the
# CRACKLES data. The data are complete, so the missing-data part masks entries
# synthetically and compares complete-case, pairwise-moment, and normal-FIML
# moment plug-ins.

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(name) name %in% args
get_val <- function(name, default, parser = as.character) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  parser(args[i + 1L])
}
if (has_flag("--help") || has_flag("-h")) {
  cat("Usage: Rscript run_experiment.R [options]\n",
      " --mask-prop P   Synthetic cell-missing fraction (default 0.25).\n",
      " --reps N        Number of masking replicates (default 25; smoke 3).\n",
      " --seed-base N   Deterministic seed base (default 26000).\n",
      " --smoke         Cheap run with 20 masking replicates.\n",
      " --help, -h      This help.\n", sep = "")
  quit("no", status = 0)
}

mask_prop <- get_val("--mask-prop", 0.25, as.numeric)
reps <- get_val("--reps", if (has_flag("--smoke")) 3L else 25L, as.integer)
seed_base <- get_val("--seed-base", 26000L, as.integer)
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

load(data_path)
d <- CRACKLES
rating_cols <- setdiff(names(d), c("patient", "UP", "LO"))
groups <- unique(sub("[0-9]+$", "", rating_cols))

ml_cov <- function(Y) {
  n <- nrow(Y)
  Z <- sweep(Y, 2L, colMeans(Y), "-")
  crossprod(Z) / n
}

make_site_layout <- function(d) {
  split_idx <- split(seq_len(nrow(d)), d$patient)
  n_by_patient <- vapply(split_idx, length, integer(1))
  if (length(unique(n_by_patient)) != 1L) {
    stop("Patients do not all have the same number of rows.", call. = FALSE)
  }
  S <- unique(n_by_patient)
  site_ord <- ave(seq_len(nrow(d)), d$patient, FUN = seq_along)
  region <- ifelse(d$UP == 1L, "UP", ifelse(d$LO == 1L, "LO", "MID"))
  region_occ <- ave(region, d$patient, region, FUN = seq_along)
  data.frame(
    row = seq_len(nrow(d)),
    patient = d$patient,
    site = site_ord,
    region = region,
    site_label = paste0(region, region_occ),
    stringsAsFactors = FALSE
  )
}

pivot_array <- function(d, cols) {
  patients <- sort(unique(d$patient))
  split_idx <- split(seq_len(nrow(d)), d$patient)
  S <- length(split_idx[[1L]])
  R <- length(cols)
  X <- array(NA_real_, dim = c(length(patients), R, S),
             dimnames = list(patient = patients, rater = cols,
                             site = paste0("site", seq_len(S))))
  for (i in seq_along(patients)) {
    rows <- split_idx[[as.character(patients[i])]]
    X[i, , ] <- t(as.matrix(d[rows, cols]))
  }
  X
}

array_to_matrix <- function(X) {
  N <- dim(X)[1L]; R <- dim(X)[2L]; S <- dim(X)[3L]
  Y <- matrix(NA_real_, nrow = N, ncol = R * S)
  cn <- character(R * S)
  k <- 0L
  for (r in seq_len(R)) {
    for (s in seq_len(S)) {
      k <- k + 1L
      Y[, k] <- X[, r, s]
      cn[k] <- paste0(dimnames(X)$rater[r], ":", dimnames(X)$site[s])
    }
  }
  colnames(Y) <- cn
  Y
}

vector_kappa <- function(mu, Sigma, R, S, W = diag(S)) {
  stopifnot(length(mu) == R * S, all(dim(Sigma) == c(R * S, R * S)))
  if (anyNA(mu) || anyNA(Sigma)) {
    return(data.frame(T = NA_real_, B = NA_real_, G = NA_real_,
                      kappa_F = NA_real_, kappa_C = NA_real_))
  }
  mu_mat <- matrix(mu, nrow = S, ncol = R)
  mu_bar <- rowMeans(mu_mat)
  T <- 0
  B <- 0
  for (r in seq_len(R)) {
    ir <- ((r - 1L) * S + 1L):(r * S)
    T <- T + sum(W * t(Sigma[ir, ir]))
    for (s in seq_len(R)) {
      is <- ((s - 1L) * S + 1L):(s * S)
      B <- B + sum(W * t(Sigma[ir, is]))
    }
  }
  G <- sum(vapply(seq_len(R), function(r) {
    delta <- mu_mat[, r] - mu_bar
    as.numeric(t(delta) %*% W %*% delta)
  }, numeric(1)))
  data.frame(
    T = T, B = B, G = G,
    kappa_F = (B - T - G) / ((R - 1) * (T + G)),
    kappa_C = (B - T) / ((R - 1) * T + R * G)
  )
}

complete_moments <- function(Y) {
  cc <- stats::complete.cases(Y)
  if (sum(cc) < 2L) {
    p <- ncol(Y)
    return(list(mu = rep(NA_real_, p),
                Sigma = matrix(NA_real_, p, p),
                n_used = sum(cc)))
  }
  list(mu = colMeans(Y[cc, , drop = FALSE]),
       Sigma = ml_cov(Y[cc, , drop = FALSE]),
       n_used = sum(cc))
}

pairwise_moments <- function(Y) {
  p <- ncol(Y)
  mu <- colMeans(Y, na.rm = TRUE)
  Sigma <- matrix(NA_real_, p, p)
  n_pair <- matrix(0L, p, p)
  for (j in seq_len(p)) {
    for (k in j:p) {
      ok <- !is.na(Y[, j]) & !is.na(Y[, k])
      n_pair[j, k] <- n_pair[k, j] <- sum(ok)
      if (sum(ok) > 0L) {
        Sigma[j, k] <- Sigma[k, j] <-
          mean((Y[ok, j] - mu[j]) * (Y[ok, k] - mu[k]))
      }
    }
  }
  Sigma[is.na(Sigma)] <- 0
  Sigma <- (Sigma + t(Sigma)) / 2
  ev <- eigen(Sigma, symmetric = TRUE)
  if (min(ev$values) < 1e-8) {
    vals <- pmax(ev$values, 1e-8)
    Sigma <- ev$vectors %*% diag(vals, nrow = length(vals)) %*% t(ev$vectors)
  }
  list(mu = mu, Sigma = Sigma, n_used = min(n_pair))
}

normal_fiml_moments <- function(Y, max_iter = 80L, tol = 1e-7) {
  N <- nrow(Y); p <- ncol(Y)
  start <- pairwise_moments(Y)
  mu <- start$mu
  Sigma <- start$Sigma + diag(1e-6, p)
  Ey_sum <- numeric(p)
  Eyy_sum <- matrix(0, p, p)
  converged <- FALSE
  iter <- 0L
  for (it in seq_len(max_iter)) {
    iter <- it
    Ey_sum[] <- 0
    Eyy_sum[,] <- 0
    for (i in seq_len(N)) {
      obs <- which(!is.na(Y[i, ]))
      mis <- which(is.na(Y[i, ]))
      yi <- mu
      V <- matrix(0, p, p)
      if (length(mis) == 0L) {
        yi <- Y[i, ]
      } else if (length(obs) == 0L) {
        V <- Sigma
      } else {
        Soo <- Sigma[obs, obs, drop = FALSE]
        Smo <- Sigma[mis, obs, drop = FALSE]
        beta <- Smo %*% solve(Soo + diag(1e-8, length(obs)))
        yi[obs] <- Y[i, obs]
        yi[mis] <- mu[mis] + beta %*% (Y[i, obs] - mu[obs])
        V[mis, mis] <- Sigma[mis, mis, drop = FALSE] -
          beta %*% Sigma[obs, mis, drop = FALSE]
      }
      Ey_sum <- Ey_sum + yi
      Eyy_sum <- Eyy_sum + V + tcrossprod(yi)
    }
    mu_new <- Ey_sum / N
    Sigma_new <- Eyy_sum / N - tcrossprod(mu_new)
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2 + diag(1e-8, p)
    delta <- max(abs(mu_new - mu), abs(Sigma_new - Sigma))
    mu <- mu_new
    Sigma <- Sigma_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }
  list(mu = mu, Sigma = Sigma, n_used = N, iter = iter, converged = converged)
}

fit_methods <- function(Y, R, S, include_fiml = TRUE) {
  methods <- list(complete = complete_moments(Y),
                  pairwise = pairwise_moments(Y))
  if (include_fiml) {
    methods$nt_fiml <- normal_fiml_moments(Y)
  }
  rows <- lapply(names(methods), function(nm) {
    m <- methods[[nm]]
    cbind(method = nm, n_used = m$n_used,
          iter = ifelse(is.null(m$iter), NA, m$iter),
          converged = ifelse(is.null(m$converged), NA, m$converged),
          vector_kappa(m$mu, m$Sigma, R, S))
  })
  do.call(rbind, rows)
}

assert_formula_checks <- function() {
  set.seed(17)
  R <- 4L; S <- 1L
  A <- matrix(rnorm(R * R), R, R)
  Sigma <- A %*% t(A)
  mu <- rnorm(R)
  got <- vector_kappa(mu, Sigma, R, S)
  T <- sum(diag(Sigma))
  B <- sum(Sigma)
  G <- sum((mu - mean(mu))^2)
  want <- (B - T - G) / ((R - 1) * (T + G))
  stopifnot(abs(got$kappa_F - want) < 1e-12)
  mu_eq <- rep(0.3, R * 3L)
  A <- matrix(rnorm((R * 3L)^2), R * 3L, R * 3L)
  Sigma <- A %*% t(A)
  got2 <- vector_kappa(mu_eq, Sigma, R, 3L)
  stopifnot(abs(got2$kappa_F - got2$kappa_C) < 1e-12)
  invisible(TRUE)
}

assert_formula_checks()

site_layout <- make_site_layout(d)
site_summary <- unique(site_layout[, c("site", "region", "site_label")])
write.csv(site_summary, file.path(results_dir, "site_layout.csv"), row.names = FALSE)

analysis_sets <- c(list(pooled = rating_cols),
                   setNames(lapply(groups, function(g) rating_cols[startsWith(rating_cols, g)]),
                            paste0("panel_", groups)))

complete_rows <- list()
masked_rows <- list()
set.seed(seed_base)
for (nm in names(analysis_sets)) {
  cols <- analysis_sets[[nm]]
  X <- pivot_array(d, cols)
  Y <- array_to_matrix(X)
  R <- dim(X)[2L]; S <- dim(X)[3L]

  full <- complete_moments(Y)
  k <- vector_kappa(full$mu, full$Sigma, R, S)
  complete_rows[[length(complete_rows) + 1L]] <- cbind(
    analysis = nm, n_patients = nrow(Y), R = R, S = S, method = "complete",
    missing_fraction = mean(is.na(Y)), n_used = full$n_used, k)

  for (rep in seq_len(reps)) {
    Ym <- Y
    mask <- matrix(stats::runif(length(Ym)) < mask_prop, nrow = nrow(Ym))
    repeat {
      Ym[mask] <- NA
      if (any(rowSums(!is.na(Ym)) == 0L) || any(colSums(!is.na(Ym)) == 0L)) {
        mask <- matrix(stats::runif(length(Ym)) < mask_prop, nrow = nrow(Ym))
      } else {
        break
      }
    }
    fits <- fit_methods(Ym, R, S, include_fiml = nm != "pooled")
    masked_rows[[length(masked_rows) + 1L]] <- cbind(
      analysis = nm, rep = rep, n_patients = nrow(Ym), R = R, S = S,
      missing_fraction = mean(is.na(Ym)), fits)
  }
}

complete_df <- do.call(rbind, complete_rows)
masked_df <- do.call(rbind, masked_rows)
for (nm in c("T", "B", "G", "kappa_F", "kappa_C", "missing_fraction")) {
  complete_df[[nm]] <- as.numeric(complete_df[[nm]])
  masked_df[[nm]] <- as.numeric(masked_df[[nm]])
}
masked_df$n_used <- as.numeric(masked_df$n_used)
masked_df$iter <- as.numeric(masked_df$iter)

write.csv(complete_df, file.path(results_dir, "complete_estimates.csv"),
          row.names = FALSE)
write.csv(masked_df, file.path(results_dir, "masked_estimates.csv"),
          row.names = FALSE)

truth <- complete_df[, c("analysis", "kappa_F", "kappa_C")]
names(truth) <- c("analysis", "truth_F", "truth_C")
masked_join <- merge(masked_df, truth, by = "analysis")
masked_join$err_F <- masked_join$kappa_F - masked_join$truth_F
masked_join$err_C <- masked_join$kappa_C - masked_join$truth_C
summarise_group <- function(g) {
  vals <- c("kappa_F", "kappa_C", "err_F", "err_C", "missing_fraction", "n_used")
  out <- g[1L, c("analysis", "method", "R", "S")]
  for (v in vals) {
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
keys <- unique(masked_join[, c("analysis", "method", "R", "S")])
summ <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  keep <- masked_join$analysis == keys$analysis[i] &
    masked_join$method == keys$method[i]
  summarise_group(masked_join[keep, ])
}))
write.csv(summ, file.path(results_dir, "masked_summary.csv"), row.names = FALSE)

metadata <- data.frame(
  key = c("mask_prop", "reps", "seed_base", "n_rows", "n_patients",
          "n_rating_columns", "R_version"),
  value = c(mask_prop, reps, seed_base, nrow(d), length(unique(d$patient)),
            length(rating_cols), paste(R.version$major, R.version$minor, sep = ".")),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(results_dir, "metadata.csv"), row.names = FALSE)

cat("Wrote:\n",
    " ", file.path(results_dir, "site_layout.csv"), "\n",
    " ", file.path(results_dir, "complete_estimates.csv"), "\n",
    " ", file.path(results_dir, "masked_estimates.csv"), "\n",
    " ", file.path(results_dir, "masked_summary.csv"), "\n",
    " ", file.path(results_dir, "metadata.csv"), "\n", sep = "")
