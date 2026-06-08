# Missingness mechanisms. Item 1 is the always-observed anchor for the MAR /
# MNAR / zero-overlap designs (keeps FIML identified). Returns X with NAs.
#
#   mcarRR   : MCAR, each cell missing w.p. RR/100 (items 2..p)
#   mcar_planned : 3-form planned missingness (structured overlap)
#   zero_overlap : split design where some item pairs are NEVER co-observed, so
#                  alpha is UNIDENTIFIED for every saturated estimator here --
#                  pairwise is undefined and saturated FIML is singular. The
#                  identification boundary: only a structural (factor) model
#                  could recover the never-observed covariances.
#   mar      : missingness on items 2..p depends on the observed anchor (item 1)
#   mnar     : missingness on item j depends on its OWN value (ignorability broken)

tune_intercept <- function(z, slope, rate) {
  f <- function(a) mean(stats::plogis(a + slope * z)) - rate
  tryCatch(stats::uniroot(f, c(-15, 15))$root, error = function(e) 0)
}

apply_mechanism <- function(X, mech, rate = 0.30) {
  n <- nrow(X); p <- ncol(X)
  if (mech == "complete") return(X)

  if (grepl("^mcar[0-9]+$", mech)) {
    rr <- as.numeric(sub("mcar", "", mech)) / 100
    tgt <- 2:p
    M <- matrix(stats::runif(n * length(tgt)) < rr, n, length(tgt))
    X[, tgt][M] <- NA_real_
    return(X)
  }

  if (mech == "mcar_planned") {
    # Anchor block = first ceil(p/4) items, always observed. Remaining items
    # split into 3 forms; each respondent randomly omits one form.
    anchor <- seq_len(max(1L, ceiling(p / 4)))
    rest <- setdiff(seq_len(p), anchor)
    if (length(rest) >= 3L) {
      forms <- split(rest, cut(seq_along(rest), 3, labels = FALSE))
      omit <- sample.int(3L, n, replace = TRUE)
      for (g in 1:3) X[omit == g, forms[[g]]] <- NA_real_
    }
    return(X)
  }

  if (mech == "zero_overlap") {
    h <- floor(p / 2)
    g1 <- setdiff(1:h, 1L)          # exclude anchor
    g2 <- setdiff((h + 1):p, 1L)
    half <- seq_len(floor(n / 2))
    X[half, g2] <- NA_real_          # first half: group-2 items missing
    X[-half, g1] <- NA_real_         # second half: group-1 items missing
    return(X)                        # pairs (a in g1, b in g2) never co-observed
  }

  if (mech == "mar") {
    z <- as.numeric(scale(X[, 1L]))
    a <- tune_intercept(z, -1.3, rate)
    for (j in 2:p) {
      pm <- stats::plogis(a + (-1.3) * z + (j - (p + 1) / 2) * 0.05)
      X[stats::runif(n) < pmin(pmax(pm, 1e-3), 0.999), j] <- NA_real_
    }
    return(X)
  }

  if (mech == "mnar") {
    for (j in 2:p) {
      z <- as.numeric(scale(X[, j]))
      a <- tune_intercept(z, -1.3, rate)
      X[stats::runif(n) < pmin(pmax(stats::plogis(a + (-1.3) * z), 1e-3), 0.999), j] <- NA_real_
    }
    return(X)
  }

  stop("unknown mechanism: ", mech)
}
