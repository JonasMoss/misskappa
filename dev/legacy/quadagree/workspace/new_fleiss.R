kappa_fun <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)

  mu <- colMeans(x)
  sigma <- cov(x) * (n - 1) / n
  a <- sum(sigma)
  mu_bar <- sum(mu) / r
  d_inv <- 1 / (sum((mu - mu_bar)^2) + sum(diag(sigma)))
  kappa <- (a * d_inv - 1) / (r - 1)

  y <- sweep(x, 2, mu, "-")
  yyt <- cbind(rowSums(y)^2, rowSums(y^2))
  y_d <- y %*% (mu - mu_bar)

  f_sigma <- c(d_inv, -a * d_inv^2)
  f_mu <- -2 * a * d_inv^2 * (mu - mu_bar)

  term_sigma <- tcrossprod(f_sigma, cov(yyt)) %*% f_sigma * (n - 1) / n
  term_mu <- tcrossprod(f_mu, sigma) %*% f_mu
  term_mu_sigma <- 4 * ((a^2 * d_inv^4) * sum(y_d * yyt[, 2])
    - (a * d_inv^3) * sum(y_d * yyt[, 1])) / n

  variance <- c((term_mu + term_mu_sigma + term_sigma) / (r - 1)^2)
  c(kappa = kappa, se = sqrt(variance) / (n-1), variance)
}

kappa_fun_2 <- \(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)

  mu <- colMeans(x)
  sigma <- cov(x) * (n - 1) / n
  a <- sum(sigma)
  mu_bar <- sum(mu) / r

  z <- cbind(y %*% (mu - mu_bar), rowSums(y)^2, rowSums(y^2))
  xi <- cov(z) * (n-1) / n

  d_fleiss <-
  d_kappa <-

  fleiss_kappa <-
  conger_kappa <-
  list(fleiss_kappa =, conger_kappa= , cov = ...)
}

kappa_fun_3 <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  r <- ncol(x)

  mu <- colMeans(x)
  sigma_mle <- cov(x) * (n - 1) / n # This is Sigma_X in math

  A0_hat <- sum(sigma_mle)
  T0_hat <- sum(diag(sigma_mle))
  mu_bar_grand <- mean(mu)
  D_mu_hat <- sum((mu - mu_bar_grand)^2) # This is sum_k (mu_k - mu_bar)^2

  # Denominators
  Den_F_hat <- (r - 1) * (D_mu_hat + T0_hat)
  Den_C_hat <- (r - 1) * T0_hat + r * D_mu_hat

  # Kappa estimates
  kappa_F_hat <- (A0_hat - (D_mu_hat + T0_hat)) / Den_F_hat
  kappa_C_hat <- (A0_hat - T0_hat) / Den_C_hat
  if(is.nan(kappa_F_hat) || is.infinite(kappa_F_hat)) kappa_F_hat <- 0 # Handle 0/0 or x/0
  if(is.nan(kappa_C_hat) || is.infinite(kappa_C_hat)) kappa_C_hat <- 0


  # Item-wise scores Z_i
  Y_centered <- sweep(x, 2, mu, "-")
  Ydm_i_scores <- Y_centered %*% (mu - mu_bar_grand) # col 1 of Z
  z1_raw_i_scores <- rowSums(Y_centered)^2          # col 2 of Z
  z2_raw_i_scores <- rowSums(Y_centered^2)          # col 3 of Z

  Z_scores_matrix <- cbind(Ydm_i_scores, z1_raw_i_scores, z2_raw_i_scores)
  Xi_hat <- cov(Z_scores_matrix) * (n - 1) / n # Cov(Z_i) with MLE scaling

  # Coefficient vectors xi_F and xi_C
  xi_F_vec <- c(2 * (1 + (r - 1) * kappa_F_hat),
                -1,
                1 + (r - 1) * kappa_F_hat)

  xi_C_vec <- c(2 * r * kappa_C_hat,
                -1,
                1 + (r - 1) * kappa_C_hat)

  # N * AVar components
  N_AVar_kappa_F <- (1 / Den_F_hat^2) * (t(xi_F_vec) %*% Xi_hat %*% xi_F_vec)
  N_AVar_kappa_C <- (1 / Den_C_hat^2) * (t(xi_C_vec) %*% Xi_hat %*% xi_C_vec)
  N_ACov_kappa_F_C <- (1 / (Den_F_hat * Den_C_hat)) * (t(xi_F_vec) %*% Xi_hat %*% xi_C_vec)

  # Asymptotic (Co)Variance matrix for (kappa_F, kappa_C)
  AVar_matrix_kappas <- matrix(c(N_AVar_kappa_F, N_ACov_kappa_F_C,
                                 N_ACov_kappa_F_C, N_AVar_kappa_C),
                               nrow = 2, ncol = 2)

  # Ensure results are single numeric values
  kappa_F_val <- as.numeric(kappa_F_hat)
  kappa_C_val <- as.numeric(kappa_C_hat)
  AVar_matrix_kappas_num <- apply(AVar_matrix_kappas, c(1,2), as.numeric)
  colnames(AVar_matrix_kappas_num) <- rownames(AVar_matrix_kappas_num) <- c("kappa_F", "kappa_C")


  list(kappa_F = kappa_F_val,
       kappa_C = kappa_C_val,
       AVar_Kappas = AVar_matrix_kappas_num)
}




zz <- cbind(y %*% (mu - mu_bar), rowSums(y)^2, rowSums(y^2) )
xi <- cov(zz) * (n-1) / n
cc <- a^2 / d^4
ff <- c(2, - d^2/a^2, 1)
ff %*% xi %*% ff * cc /  (r - 1)^2

gg <- c(-2*a/d^2, 1/d, -a/d^2)
gg %*% xi %*% gg /  (r - 1)^2

gg <- c(2*a/d, -1, a/d)
gg %*% xi %*% gg /  (r - 1)^2 / d^2

gg <- c(2*a, -d, a)
gg %*% xi %*% gg /  (r - 1)^2 / d^4

gg <- c(2, -d / a, 1)
gg %*% xi %*% gg /  (r - 1)^2 * a^2 / d^4

xi_ <- a/d
gg <- c(2*xi_, -1, xi_)
gg %*% xi %*% gg /  (r - 1)^2  / d^2


kappa_c <- conger(x)

gg <- c(2*r*kappa_c, -1, 1 + (r-1)*kappa_c)
gg %*% xi %*% gg




x <- as.matrix(x)
n <- nrow(x)
r <- ncol(x)

mu <- colMeans(x)
sigma <- cov(x) * (n - 1) / n
a <- sum(sigma)
mu_bar <- sum(mu) / r
d_inv <- 1 / (sum((mu - mu_bar)^2) + sum(diag(sigma)))
kappa <- (a * d_inv - 1) / (r - 1)

y <- sweep(x, 2, mu, "-")
yyt <- cbind(rowSums(y)^2, rowSums(y^2))
y_d <- y %*% (mu - mu_bar)

f_sigma <- c(d_inv, -a * d_inv^2)
f_mu <- -2 * a * d_inv^2 * (mu - mu_bar)

term_sigma <- tcrossprod(f_sigma, cov(yyt)) %*% f_sigma * (n - 1) / n
term_mu <- tcrossprod(f_mu, sigma) %*% f_mu
term_mu_sigma <- 4 * ((a^2 * d_inv^4) * sum(y_d * yyt[, 2]) - (a * d_inv^3) * sum(y_d * yyt[, 1])) / n

avar <- c((term_mu + term_mu_sigma + term_sigma) / (r - 1)^2)









# yyt <- sapply(seq(n), \(i) {
#   yyt <- tcrossprod(y[i, ])
#   c(sum(yyt), sum(diag(yyt)))
# })
#
# yyt2 <- cbind(apply(y, 1,  \(y) sum(tcrossprod(y))), rowSums(y^2))
# f <- c(1, -a/d)
# tcrossprod(f %*% cov(t(yyt)), f) / d^2 * (n-1) / n
#
# z1_i_vec <- numeric(n)
# z2_i_vec <- numeric(n)
#
# for (i_item in 1:n) {
#   item <- tcrossprod(y[i_item, ]) # R x R matrix for item i: Yi %*% t(Yi)
#   z1_i_vec[i_item] <- sum(item)
#   z2_i_vec[i_item] <- sum(diag(item))
# }
#
# var_z1 <- var(z1_i_vec) * (n-1)/n
# var_z2 <- var(z2_i_vec) * (n-1)/n
# cov_z1_z2 <- cov(z1_i_vec, z2_i_vec) * (n-1)/n
#
# # Then f_sigma^T Gamma f_sigma term is:
# term2_var_f_star_new <- Delta_R_hat^(-2) * var_z1 +
#   a^2 * Delta_R_hat^(-4) * var_z2 -
#   2 * a * Delta_R_hat^(-3) * cov_z1_z2
#
# term2_var_f_star_new <- d^(-2) * var_z1 +
#   a^2 * d^(-4) * var_z2 -
#   2 * a * d^(-3) * cov_z1_z2
#
#
# cov(t(zz)) * (n-1) / n
#
# f <- c(1, -a/d)
# tcrossprod(f %*% cov(t(zz)) * (n-1) / n, f) / d^2
#
