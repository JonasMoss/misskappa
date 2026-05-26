conger <- function(ratings, weights = "unweighted", categ.labels = NULL,
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
      weights.mat <- irrCAC::quadratic.weights(categ)
    } else if (weights == "ordinal") {
      weights.mat <- irrCAC::ordinal.weights(categ)
    } else if (weights == "linear") {
      weights.mat <- irrCAC::linear.weights(categ)
    } else if (weights == "radical") {
      weights.mat <- irrCAC::radical.weights(categ)
    } else if (weights == "ratio") {
      weights.mat <- irrCAC::ratio.weights(categ)
    } else if (weights == "circular") {
      weights.mat <- irrCAC::circular.weights(categ)
    } else if (weights == "bipolar") {
      weights.mat <- irrCAC::bipolar.weights(categ)
    } else {
      weights.mat <- irrCAC::identity.weights(categ)
    }
  } else {
    w.name <- "Custom Weights"
    weights.mat <- weights
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
  classif.mat <- matrix(0, nrow = r, ncol = q)
  for (k in 1:q) {
    with.mis <- (t(ratings.mat) == categ[k])
    without.mis <- replace(with.mis, is.na(with.mis), FALSE)
    classif.mat[, k] <- without.mis %*% rep(1, n)
  }
  ri.vec <- agree.mat %*% rep(1, q)
  sum.q <- (agree.mat * (agree.mat.w - 1)) %*% rep(1, q)
  n2more <- sum(ri.vec >= 2)
  pa <- sum(sum.q[ri.vec >= 2] / ((ri.vec * (ri.vec - 1))[ri.vec >=
    2])) / n2more
  ng.vec <- classif.mat %*% rep(1, q)
  pgk.mat <- classif.mat / (ng.vec %*% rep(1, q))
  p.mean.k <- (t(pgk.mat) %*% rep(1, r)) / r
  s2kl.mat <- (t(pgk.mat) %*% pgk.mat - r * p.mean.k %*% t(p.mean.k)) / (r -
    1)
  pe <- sum(weights.mat * (p.mean.k %*% t(p.mean.k) - s2kl.mat / r))
  conger.kappa <- (pa - pe) / (1 - pe)
  conger.kappa.est <- conger.kappa
  bkl.mat <- (weights.mat + t(weights.mat)) / 2
  pe.ivec1 <- r * (agree.mat %*% t(t(p.mean.k) %*% bkl.mat))
  pe.ivec2 <- rep(0, n)
  lamda.ig.mat <- matrix(0, n, r)
  if (is.numeric(ratings.mat)) {
    epsi.ig.mat <- 1 - is.na(ratings.mat)
    epsi.ig.mat <- replace(
      epsi.ig.mat, is.na(epsi.ig.mat),
      FALSE
    )
  } else {
    epsi.ig.mat <- 1 - (ratings.mat == "")
    epsi.ig.mat <- replace(
      epsi.ig.mat, is.na(epsi.ig.mat),
      FALSE
    )
  }
  for (k in 1:q) {
    lamda.ig.kmat <- matrix(0, n, r)
    for (l in 1:q) {
      delta.ig.mat <- (ratings.mat == categ[l])
      delta.ig.mat <- replace(
        delta.ig.mat, is.na(delta.ig.mat),
        FALSE
      )
      lamda.ig.kmat <- lamda.ig.kmat + weights.mat[k, l] *
        (delta.ig.mat - (epsi.ig.mat - rep(1, n) %*%
          t(ng.vec / n)) * (rep(1, n) %*% t(pgk.mat[, l])))
    }
    lamda.ig.kmat <- lamda.ig.kmat * (rep(1, n) %*% t(n / ng.vec))
    lamda.ig.mat <- lamda.ig.mat + lamda.ig.kmat * (r * mean(pgk.mat[
      ,
      k
    ]) - rep(1, n) %*% t(pgk.mat[, k]))
  }
  pe.ivec <- (lamda.ig.mat %*% rep(1, r)) / (r * (r - 1))
  den.ivec <- ri.vec * (ri.vec - 1)
  den.ivec <- den.ivec - (den.ivec == 0)
  pa.ivec <- sum.q / den.ivec
  pe.r2 <- pe * (ri.vec >= 2)
  conger.ivec <- (n / n2more) * (pa.ivec - pe.r2) / (1 - pe)
  conger.ivec.x <- conger.ivec - 2 * (1 - conger.kappa) * (pe.ivec -
    pe) / (1 - pe)
  var.conger <- NA
  stderr <- NA
  stderr.est <- NA
  p.value <- NA
  lcb <- NA
  ucb <- NA
  if (n >= 2) {
    var.conger <- ((1 - f) / (n * n)) * sum((conger.ivec.x -
      conger.kappa)^2)
    stderr <- sqrt(var.conger)
    p.value <- 2 * (1 - pt(abs(conger.kappa / stderr), n -
      1))
    lcb <- conger.kappa - stderr * qt(
      1 - (1 - conflev) / 2,
      n - 1
    )
    ucb <- min(1, conger.kappa + stderr * qt(
      1 - (1 - conflev) / 2,
      n - 1
    ))
  }

  c(kappa = conger.kappa, stderr = stderr)
}

fleiss <- function(ratings, weights = "unweighted", categ = NULL, conflev = 0.95,
                   N = Inf) {
  agree.mat <- as.matrix(ratings)
  n <- nrow(agree.mat)
  q <- ncol(agree.mat)
  f <- n / N
  if (is.null(categ)) {
    categ <- 1:q
  } else {
    q2 <- length(categ)
    if (!is.numeric(categ)) {
      categ <- 1:q2
    }
    if (q2 > q) {
      colna1 <- colnames(agree.mat)
      agree.mat <- cbind(agree.mat, matrix(0, n, (q2 -
        q)))
      colna2 <- sapply(1:(q2 - q), function(x) {
        paste0(
          "v",
          x
        )
      })
      colnames(agree.mat) <- c(colna1, colna2)
      q <- q2
    }
  }
  if (is.character(weights)) {
    w.name <- weights
    if (weights == "quadratic") {
      weights.mat <- irrCAC::quadratic.weights(categ)
    } else if (weights == "ordinal") {
      weights.mat <- irrCAC::ordinal.weights(categ)
    } else if (weights == "linear") {
      weights.mat <- irrCAC::linear.weights(categ)
    } else if (weights == "radical") {
      weights.mat <- irrCAC::radical.weights(categ)
    } else if (weights == "ratio") {
      weights.mat <- irrCAC::ratio.weights(categ)
    } else if (weights == "circular") {
      weights.mat <- irrCAC::circular.weights(categ)
    } else if (weights == "bipolar") {
      weights.mat <- irrCAC::bipolar.weights(categ)
    } else {
      weights.mat <- irrCAC::identity.weights(categ)
    }
  } else {
    w.name <- "Custom Weights"
    weights.mat <- as.matrix(weights)
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
  var.fleiss <- ((1 - f) / (n * n)) * sum((kappa.ivec.x -
    fleiss.kappa)^2)
  stderr <- sqrt(var.fleiss)
  p.value <- 2 * (1 - pt(fleiss.kappa / stderr, n - 1))
  lcb <- fleiss.kappa - stderr * qt(1 - (1 - conflev) / 2, n -
    1)
  ucb <- min(1, fleiss.kappa + stderr * qt(
    1 - (1 - conflev) / 2,
    n - 1
  ))
  conf.int <- paste0(
    "(", round(lcb, 3), ",", round(ucb, 3),
    ")"
  )
  coeff <- fleiss.kappa
  coeff.name <- "Fleiss' Kappa"
  c(kappa = fleiss.kappa, stderr = stderr)
}

bp <- function(ratings, weights = "unweighted", categ = NULL, conflev = 0.95,
               N = Inf) {
  agree.mat <- as.matrix(ratings)
  n <- nrow(agree.mat)
  q <- ncol(agree.mat)
  f <- n / N
  if (is.null(categ)) {
    categ <- 1:q
  } else {
    q2 <- length(categ)
    if (!is.numeric(categ)) {
      categ <- 1:q2
    }
    if (q2 > q) {
      colna1 <- colnames(agree.mat)
      agree.mat <- cbind(agree.mat, matrix(0, n, (q2 -
        q)))
      colna2 <- sapply(1:(q2 - q), function(x) {
        paste0(
          "v",
          x
        )
      })
      colnames(agree.mat) <- c(colna1, colna2)
      q <- q2
    }
  }
  if (is.character(weights)) {
    w.name <- weights
    if (weights == "quadratic") {
      weights.mat <- irrCAC::quadratic.weights(categ)
    } else if (weights == "ordinal") {
      weights.mat <- irrCAC::ordinal.weights(categ)
    } else if (weights == "linear") {
      weights.mat <- irrCAC::linear.weights(categ)
    } else if (weights == "radical") {
      weights.mat <- irrCAC::radical.weights(categ)
    } else if (weights == "ratio") {
      weights.mat <- irrCAC::ratio.weights(categ)
    } else if (weights == "circular") {
      weights.mat <- irrCAC::circular.weights(categ)
    } else if (weights == "bipolar") {
      weights.mat <- irrCAC::bipolar.weights(categ)
    } else {
      weights.mat <- irrCAC::identity.weights(categ)
    }
  } else {
    w.name <- "Custom Weights"
    weights.mat <- as.matrix(weights)
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
  pe <- sum(weights.mat) / (q^2)
  bp.coeff <- (pa - pe) / (1 - pe)
  den.ivec <- ri.vec * (ri.vec - 1)
  den.ivec <- den.ivec - (den.ivec == 0)
  pa.ivec <- sum.q / den.ivec
  pe.r2 <- pe * (ri.vec >= 2)
  bp.ivec <- (n / n2more) * (pa.ivec - pe.r2) / (1 - pe)
  var.bp <- ((1 - f) / (n * n)) * sum((bp.ivec - bp.coeff)^2)
  stderr <- sqrt(var.bp)
  p.value <- 2 * (1 - pt(bp.coeff / stderr, n - 1))
  lcb <- bp.coeff - stderr * qt(1 - (1 - conflev) / 2, n - 1)
  ucb <- min(1, bp.coeff + stderr * qt(
    1 - (1 - conflev) / 2,
    n - 1
  ))
  conf.int <- paste0(
    "(", round(lcb, 3), ",", round(ucb, 3),
    ")"
  )
  coeff <- bp.coeff
  coeff.name <- "Brennan-Prediger"
  df.out <- data.frame(
    coeff.name, coeff, stderr, conf.int,
    p.value, pa, pe
  )

  c(kappa = bp.coeff, stderr = stderr)
}
