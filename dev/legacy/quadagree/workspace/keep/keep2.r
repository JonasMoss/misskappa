fleiss.kappa.unrounded_from_irrCAC <- function(ratings, weights = "unweighted", categ.labels = NULL,
                                               conflev = 0.95, N = Inf) {
  ratings.mat <- as.matrix(ratings)
  if (is.character(ratings.mat)) {
    ratings.mat.upper <- toupper(ratings.mat)
    ratings.mat.trimmed <- trimws(ratings.mat.upper)
    ratings.mat.trimmed[ratings.mat.trimmed == ""] <- NA_character_
    ratings.mat <- ratings.mat.trimmed
  }
  n <- nrow(ratings.mat)
  r <- ncol(ratings.mat)
  f <- n / N
  if (is.null(categ.labels)) {
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  } else {
    categ <- toupper(categ.labels)
  }
  q <- length(categ)

  # Helper functions from irrCAC (or their equivalents if not exported)
  # Assuming these exist or are simple enough to replicate if needed.
  # For this example, I'll assume irrCAC's internal weight functions would be called.
  # We need to be careful if they aren't exported.
  # Let's stub them for now if direct calls like irrCAC:::quadratic.weights aren't viable
  # or if we want this to be truly standalone from irrCAC's non-exported parts.

  # Replicating irrCAC's weight matrix logic broadly
  if (is.character(weights)) {
    w.name <- weights
    if (weights == "quadratic") {
      # Simplified quadratic weights logic - irrCAC's is more general
      if (is.numeric(categ)) {
        dist.mat <- abs(outer(categ, categ, "-"))
        max.dist <- max(dist.mat)
        weights.mat <- 1 - (dist.mat / max.dist)^2
      } else { # Fallback for non-numeric, assumes identity if not handled
        weights.mat <- diag(q)
      }
    } else { # Default to identity weights if not quadratic
      weights.mat <- diag(q)
      if (weights != "unweighted") {
        w.name <- "Custom (Identity Fallback)"
      } else {
        w.name <- "unweighted"
      }
    }
  } else {
    w.name <- "Custom Weights"
    weights.mat <- as.matrix(weights)
  }
  if (weights == "unweighted") weights.mat <- diag(q) # Ensure identity for unweighted

  agree.mat <- matrix(0, nrow = n, ncol = q)
  for (k in 1:q) {
    categ.is.k <- (ratings.mat == categ[k])
    agree.mat[, k] <- rowSums(replace(categ.is.k, is.na(categ.is.k), FALSE))
  }

  agree.mat.w <- t(weights.mat %*% t(agree.mat))
  ri.vec <- rowSums(agree.mat)

  valid_raters_per_item <- ri.vec * (ri.vec - 1)
  items_with_sufficient_raters <- valid_raters_per_item > 0

  n2more <- sum(items_with_sufficient_raters)

  if (n2more == 0) { # Handle case with no items having >=2 ratings
    pa <- 0
  } else {
    sum_q_val <- rowSums(agree.mat * (agree.mat.w - agree.mat %*% (1 - weights.mat)))
    pa_num <- sum_q_val[items_with_sufficient_raters]
    pa_den <- valid_raters_per_item[items_with_sufficient_raters]
    pa <- sum(pa_num / pa_den) / n2more
  }

  pi.vec_num <- colSums(agree.mat, na.rm = TRUE)
  pi.vec_den <- sum(ri.vec, na.rm = TRUE)
  if (pi.vec_den == 0) pi.vec <- rep(0, q) else pi.vec <- pi.vec_num / pi.vec_den


  pe <- sum(weights.mat * (pi.vec %*% t(pi.vec)))
  fleiss.kappa <- (pa - pe) / (1 - pe)

  # --- Start of Variance Calculation (largely from irrCAC) ---
  den.ivec <- ri.vec * (ri.vec - 1)
  den.ivec.safe <- den.ivec
  den.ivec.safe[den.ivec == 0] <- 1 # Avoid division by zero, result will be NA or handled

  sum_q_for_pa_ivec <- rowSums(agree.mat * (agree.mat.w - agree.mat %*% (1 - weights.mat)))

  pa.ivec <- sum_q_for_pa_ivec / den.ivec.safe
  pa.ivec[den.ivec == 0] <- 0 # Or NA, depending on how irrCAC handles it internally

  pe.r2 <- pe * (ri.vec >= 2) # This is a conditional pe for items with >=2 ratings

  # kappa.ivec needs careful construction if there are items with < 2 ratings.
  # irrCAC uses n/n2more scaling factor.
  kappa.ivec_num <- (pa.ivec - pe.r2)
  kappa.ivec_den <- (1 - pe)
  if (kappa.ivec_den == 0) kappa.ivec_den <- 1e-9 # Avoid division by zero if 1-pe is 0

  # Only scale kappa.ivec for items that contribute to n2more
  kappa.ivec <- rep(0, n)
  kappa.ivec[items_with_sufficient_raters] <- (n / n2more) *
    kappa.ivec_num[items_with_sufficient_raters] /
    kappa.ivec_den


  pi.vec.wk. <- weights.mat %*% pi.vec
  pi.vec.w.k <- t(weights.mat) %*% pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k) / 2

  pe.ivec_num <- agree.mat %*% pi.vec.w
  pe.ivec_den <- ri.vec
  pe.ivec_den_safe <- pe.ivec_den
  pe.ivec_den_safe[pe.ivec_den_safe == 0] <- 1 # Avoid division by zero
  pe.ivec <- pe.ivec_num / pe.ivec_den_safe
  pe.ivec[pe.ivec_den == 0] <- pe # if no raters for item, pe.ivec might be pe? irrCAC is complex here

  kappa.ivec.x_num <- (pe.ivec - pe)
  kappa.ivec.x_den <- (1 - pe)
  if (kappa.ivec.x_den == 0) kappa.ivec.x_den <- 1e-9

  kappa.ivec.x <- kappa.ivec - 2 * (1 - fleiss.kappa) * kappa.ivec.x_num / kappa.ivec.x_den

  var.fleiss <- NA
  stderr <- NA

  if (n >= 2 && n2more > 0 && (n - 1) > 0) { # Added n2more > 0 and n-1 > 0 check
    var.fleiss <- ((1 - f) / (n * (n - 1))) * sum((kappa.ivec.x - fleiss.kappa)^2, na.rm = TRUE)
    # Ensure var.fleiss is not negative due to precision issues, though unlikely here
    if (!is.na(var.fleiss) && var.fleiss < 0) var.fleiss <- 0
    stderr <- sqrt(var.fleiss)
  }
  # --- End of Variance Calculation ---

  # Return unrounded values
  return(list(
    kappa_unrounded = fleiss.kappa,
    stderr_unrounded = stderr,
    pa = pa,
    pe = pe,
    n_items = n,
    n_raters = r,
    n_categories = q,
    n_items_ge2_ratings = n2more
  ))
}

# Helper for trimws if not on new enough R
if (!exists("trimws")) {
  trimws <- function(x, which = c("both", "left", "right")) {
    which <- match.arg(which)
    if (which == "both") {
      return(gsub("^[ \t\r\n]+|[ \t\r\n]+$", "", x))
    }
    if (which == "left") {
      return(gsub("^[ \t\t\r\n]+", "", x))
    }
    if (which == "right") {
      return(gsub("[ \t\r\n]+$", "", x))
    }
  }
}

ratings <- dat.zapf2016
weights <- "quadratic"
categ.labels <- NULL
conflev <- 0.95
N <- Inf



irrnew <- function(ratings, weights = "quadratic", categ.labels = NULL,
                   conflev = 0.95, N = Inf) {
  ratings.mat <- as.matrix(ratings)
  if (is.character(ratings.mat)) {
    ratings.mat <- trim(toupper(ratings.mat))
    ratings.mat[ratings.mat == ""] <- NA_character_
  }
  n <- nrow(ratings.mat)
  r <- ncol(ratings.mat)
  f <- n / N
  if (is.null(categ.labels)) {
    categ.init <- unique(na.omit(as.vector(ratings.mat)))
    categ <- sort(categ.init)
  } else {
    categ <- toupper(categ.labels)
  }
  q <- length(categ)
  if (is.character(weights)) {
    w.name <- weights
    if (weights == "quadratic") {
      weights.mat <- irrCAC:::quadratic.weights(categ)
    } else if (weights == "ordinal") {
      weights.mat <- ordinal.weights(categ)
    } else if (weights == "linear") {
      weights.mat <- linear.weights(categ)
    } else if (weights == "radical") {
      weights.mat <- radical.weights(categ)
    } else if (weights == "ratio") {
      weights.mat <- ratio.weights(categ)
    } else if (weights == "circular") {
      weights.mat <- circular.weights(categ)
    } else if (weights == "bipolar") {
      weights.mat <- bipolar.weights(categ)
    } else {
      weights.mat <- identity.weights(categ)
    }
  } else {
    w.name <- "Custom Weights"
    weights.mat <- as.matrix(weights)
  }
  agree.mat <- matrix(0, nrow = n, ncol = q)
  for (k in 1:q) {
    categ.is.k <- (ratings.mat == categ[k])
    agree.mat[, k] <- (replace(
      categ.is.k, is.na(categ.is.k),
      FALSE
    )) %*% rep(1, r)
  }

  agree.mat.w <- t(weights.mat %*% t(agree.mat))
  ri.vec <- agree.mat %*% rep(1, q)
  sum.q <- (agree.mat * (agree.mat.w - 1)) %*% rep(1, q)
  n2more <- sum(ri.vec >= 2)
  pa <- sum(sum.q[ri.vec >= 2] / ((ri.vec * (ri.vec - 1))[ri.vec >=
    2])) / n2more
  pi.vec <- t(t(rep(1 / n, n)) %*% (agree.mat / (ri.vec %*% t(rep(
    1,
    q
  )))))
  pe <- sum(weights.mat * (pi.vec %*% t(pi.vec)))
  fleiss.kappa <- (pa - pe) / (1 - pe)
  fleiss.kappa.est <- fleiss.kappa
  den.ivec <- ri.vec * (ri.vec - 1)
  den.ivec <- den.ivec - (den.ivec == 0)
  pa.ivec <- sum.q / den.ivec
  pe.r2 <- pe * (ri.vec >= 2)
  kappa.ivec <- (n / n2more) * (pa.ivec - pe.r2) / (1 - pe)
  pi.vec.wk. <- weights.mat %*% pi.vec
  pi.vec.w.k <- t(weights.mat) %*% pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k) / 2
  pe.ivec <- (agree.mat %*% pi.vec.w) / ri.vec
  kappa.ivec.x <- kappa.ivec - 2 * (1 - fleiss.kappa) * (pe.ivec -
    pe) / (1 - pe)
  var.fleiss <- NA
  stderr <- NA
  stderr.est <- NA
  p.value <- NA
  lcb <- NA
  ucb <- NA
  var.fleiss <- (1 / (n * (n - 1))) * sum((kappa.ivec.x - fleiss.kappa)^2)
  sqrt(var.fleiss)
}

irrnew(dat.zapf2016, weights = "quadratic") * sqrt(n-1)
sqrt(get_fk_var_robust_full(dat.zapf2016) * nrow(dat.zapf2016))
irrCAC::fleiss.kappa.raw(dat.zapf2016, weights = "quadratic")
fleiss.kappa.unrounded_from_irrCAC(dat.zapf2016, weights = "quadratic")
