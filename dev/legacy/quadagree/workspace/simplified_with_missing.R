get_biased_pairwise_cov <- function(x) {
  r <- ncol(x)
  n <- nrow(x)
  cov_mat <- matrix(0, r, r)
  for (j in 1:r) {
    for (k in j:r) {
      complete_cases <- !is.na(x[, j]) & !is.na(x[, k])
      n_jk <- sum(complete_cases)
      if (n_jk > 0) {
        x_j <- x[complete_cases, j]
        x_k <- x[complete_cases, k]
        # Use cov.wt from {stats} which allows specifying the denominator
        # The default weights are 1, so this calculates sum of products / n_jk
        # This is equivalent to cross-product of centered variables / n_jk
        cov_mat[j, k] <- cov.wt(cbind(x_j, x_k), method = "ML")$cov[1, 2]
      }
    }
  }
  # Symmetrize
  cov_mat[lower.tri(cov_mat)] <- t(cov_mat)[lower.tri(cov_mat)]
  return(cov_mat)
}

.get_moments_and_probs <- function(x) {
  n <- nrow(x)
  r <- ncol(x)
  m_vec <- !is.na(x)

  mu_hat <- colMeans(x, na.rm = TRUE)
  y_hat <- sweep(x, 2, mu_hat, "-")

  p1 <- colMeans(m_vec)
  p2 <- crossprod(m_vec) / n

  mu3 <- array(0, dim = c(r, r, r))
  mu4 <- array(0, dim = c(r, r, r, r))
  p3 <- array(0, dim = c(r, r, r))
  p4 <- array(0, dim = c(r, r, r, r))

  for (i in 1:r) {
    for (j in 1:r) {
      for (k in 1:r) {
        mu3[i, j, k] <- mean(y_hat[, i] * y_hat[, j] * y_hat[, k], na.rm = TRUE)
        p3[i, j, k] <- mean(m_vec[, i] * m_vec[, j] * m_vec[, k])
        for (l in 1:r) {
          mu4[i, j, k, l] <- mean(y_hat[, i] * y_hat[, j] * y_hat[, k] * y_hat[, l], na.rm = TRUE)
          p4[i, j, k, l] <- mean(m_vec[, i] * m_vec[, j] * m_vec[, k] * m_vec[, l])
        }
      }
    }
  }

  list(
    r = r,
    mu_hat = mu_hat,
    sigma_hat = get_biased_pairwise_cov(x),
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    mu3 = mu3,
    mu4 = mu4
  )
}

.get_psi_matrix <- function(d) {
  r <- d$r
  mu_hat <- d$mu_hat
  sigma_hat <- d$sigma_hat
  p1 <- d$p1
  p2 <- d$p2
  p3 <- d$p3
  p4 <- d$p4
  mu3 <- d$mu3
  mu4 <- d$mu4

  Psi <- matrix(NA, 3, 3)

  get_cov_of_sums <- function(v1, v2, type) {
    total_cov <- 0
    if (type == "sigma-sigma") {
      for (i in 1:r) {
        for (j in 1:r) {
          for (k in 1:r) {
            for (l in 1:r) {
              gamma <- mu4[i, j, k, l] - sigma_hat[i, j] * sigma_hat[k, l]
              pi_gamma <- p4[i, j, k, l] / (p2[i, j] * p2[k, l])
              total_cov <- total_cov + v1[i, j] * v2[k, l] * pi_gamma * gamma
            }
          }
        }
      }
    } else if (type == "mu-mu") {
      for (i in 1:r) {
        for (j in 1:r) {
          omega_ij <- (p2[i, j] / (p1[i] * p1[j])) * sigma_hat[i, j]
          total_cov <- total_cov + v1[i] * v2[j] * omega_ij
        }
      }
    } else if (type == "mu-sigma") {
      for (i in 1:r) {
        for (j in 1:r) {
          for (k in 1:r) {
            omega_ijk <- (p3[i, j, k] / (p1[i] * p2[j, k])) * mu3[i, j, k]
            total_cov <- total_cov + v1[i] * v2[j, k] * omega_ijk
          }
        }
      }
    }
    return(total_cov)
  }

  v_total <- matrix(1, r, r)
  v_trace <- diag(1, r, r)
  v_dmu <- 2 * (mu_hat - mean(mu_hat))

  Psi[1, 1] <- get_cov_of_sums(v_total, v_total, "sigma-sigma")
  Psi[2, 2] <- get_cov_of_sums(v_trace, v_trace, "sigma-sigma")
  Psi[3, 3] <- get_cov_of_sums(v_dmu, v_dmu, "mu-mu")
  Psi[1, 2] <- Psi[2, 1] <- get_cov_of_sums(v_total, v_trace, "sigma-sigma")
  Psi[1, 3] <- Psi[3, 1] <- get_cov_of_sums(v_dmu, v_total, "mu-sigma")
  Psi[2, 3] <- Psi[3, 2] <- get_cov_of_sums(v_dmu, v_trace, "mu-sigma")

  return(Psi)
}

.get_gradient_matrix <- function(d) {
  r <- d$r
  sigma_hat <- d$sigma_hat
  mu_hat <- d$mu_hat

  theta1_h <- sum(sigma_hat)
  theta2_h <- sum(diag(sigma_hat))
  theta3_h <- sum((mu_hat - mean(mu_hat))^2)

  NF <- theta1_h - theta2_h - theta3_h
  DF <- (r - 1) * theta2_h + (r - 1) * theta3_h
  grad_F <- c(1 / DF, (-DF - NF * (r - 1)) / DF^2, (-DF - NF * (r - 1)) / DF^2)

  NC <- theta1_h - theta2_h
  DC <- (r - 1) * theta2_h + r * theta3_h
  grad_C <- c(1 / DC, (-DC - NC * (r - 1)) / DC^2, (-NC * r) / DC^2)

  return(cbind(grad_F, grad_C))
}

kappa_acov_aggregate <- function(x) {
  x <- as.matrix(x)
  d <- .get_moments_and_probs(x)

  Psi <- .get_psi_matrix(d)
  G <- .get_gradient_matrix(d)

  scaled_acov <- t(G) %*% Psi %*% G

  tr_s <- sum(diag(d$sigma_hat))
  d_mu <- sum((d$mu_hat - mean(d$mu_hat))^2)

  NF <- sum(d$sigma_hat) - tr_s - d_mu
  DF <- (d$r - 1) * tr_s + (d$r - 1) * d_mu
  NC <- sum(d$sigma_hat) - tr_s
  DC <- (d$r - 1) * tr_s + d$r * d_mu

  list(
    fleiss = NF / DF,
    conger = NC / DC,
    acov = scaled_acov / nrow(x),
    scaled_acov
  )
}


# --- Verification Run ---
# We run this new function against our trusted formal one.
set.seed(135)
n_sim <- 50
r_sim <- 4
missing_prob <- 0.2
true_mu <- c(1, 1.5, 1.2, 2)
true_sigma <- matrix(0.7, r_sim, r_sim)
diag(true_sigma) <- 1
L <- chol(true_sigma)
complete_data <- matrix(rnorm(n_sim * r_sim), n_sim, r_sim) %*% L
complete_data <- sweep(complete_data, 2, true_mu, "+")
x_missing <- complete_data
x_missing[sample(1:(n_sim * r_sim), size = n_sim * r_sim * missing_prob)] <- NA

## Test
kappa_acov_aggregate(x_missing)$acov[2,2]
#congerci_rcpp(x_missing)$sd^2
kappa_raw_rcpp(x_missing)$scaled_acov[2, 2] / kappa_raw_rcpp(x_missing)$n_eff

kappa_raw_rcpp(x_missing)
kappa_raw(x_missing)
KappasRawCpp(complete_data, weight="quadratic")


irrCAC::fleiss.kappa.raw(complete_data, weight="quadratic")$est
kappa_raw(complete_data)
#kappa_raw_rcpp(x_missing)$scaled_acov / kappa_raw_rcpp(x_missing)$n_eff

#kappa_raw_rcpp(x_missing)


