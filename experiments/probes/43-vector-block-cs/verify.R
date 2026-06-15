# Probe 43 — vector-valued block-compound-symmetry (block-CS) covariance for the
# quadratic Conger/Fleiss/BP coefficients. Before the C++ surgery, lock the two
# mathematical claims the implementation rests on:
#
# (1) PROJECTION = MLE. For block-CS Sigma = J_R (x) Psi + I_R (x) E (Psi, E both
#     d x d), the constrained Gaussian MLE given the unconstrained covariance S is
#     block averaging: Psi_hat = mean off-diagonal d x d block of S,
#     E_hat = mean diagonal block - Psi_hat. (Same fixed rater-eigenbasis argument
#     as the scalar case; d=1 reduces to project_variance_components.)
#
# (2) COEFFICIENT CONSISTENCY. The vector quadratic Conger/Fleiss/BP computed from
#     (mu, block-CS Sigma) via the J/I/(I-J/R) (x) W contractions are consistent
#     for the population coefficients under MCAR.

set.seed(20260615)

# ---- block-CS building blocks ------------------------------------------------
make_sigma <- function(R, Psi, E) {                 # J_R (x) Psi + I_R (x) E
  d <- nrow(Psi)
  S <- kronecker(matrix(1, R, R), Psi) + kronecker(diag(R), E)
  S
}
block <- function(S, r, s, d) S[((r - 1) * d + 1):(r * d), ((s - 1) * d + 1):(s * d)]

project_block_cs <- function(S, R, d) {             # block averaging
  Dbar <- matrix(0, d, d); Obar <- matrix(0, d, d)
  for (r in 1:R) Dbar <- Dbar + block(S, r, r, d)
  Dbar <- Dbar / R
  for (r in 1:R) for (s in 1:R) if (r != s) Obar <- Obar + block(S, r, s, d)
  Obar <- Obar / (R * (R - 1))
  Psi <- (Obar + t(Obar)) / 2
  E <- (Dbar - Psi); E <- (E + t(E)) / 2
  list(Psi = Psi, E = E)
}

# vector quadratic coefficients from (mu stacked Rd, Sigma Rd x Rd, W d x d)
vec_coef <- function(mu, Sigma, R, d, W, bp_chance) {
  Bd <- kronecker(matrix(1, R, R), W)
  Td <- kronecker(diag(R), W)
  Gd <- kronecker(diag(R) - matrix(1 / R, R, R), W)
  tB <- sum(Bd * Sigma); tT <- sum(Td * Sigma); tG <- as.numeric(t(mu) %*% Gd %*% mu)
  conger <- (tB - tT) / ((R - 1) * tT + R * tG)
  fleiss <- (tB - tT - tG) / ((R - 1) * (tT + tG))
  d_obs <- (2 / (R - 1)) * (tT + tG - tB / R)
  bp <- 1 - d_obs / bp_chance
  c(Conger = conger, Fleiss = fleiss, `Brennan-Prediger` = bp)
}

# ---- (1) projection == MLE on complete data ----------------------------------
cat("=== (1) block-averaging projection vs numerical Gaussian MLE (complete data) ===\n")
R <- 4; d <- 2
Psi_true <- matrix(c(1, .4, .4, .8), 2, 2)
E_true   <- matrix(c(.6, .1, .1, .5), 2, 2)
mu_true  <- c(0, .2,  .5, .1,  1, -.3,  .7, .4)        # R*d = 8 free means
Sig_true <- make_sigma(R, Psi_true, E_true)

n <- 8000
L <- chol(Sig_true)
Y <- matrix(mu_true, n, R * d, byrow = TRUE) + matrix(rnorm(n * R * d), n, R * d) %*% L
mu_hat <- colMeans(Y)
S <- crossprod(scale(Y, center = TRUE, scale = FALSE)) / n
pj <- project_block_cs(S, R, d)

# numerical MLE of the Gaussian over (Psi, E) at mu = mu_hat (S is the MLE cov)
negll <- function(par) {
  Lp <- matrix(0, d, d); Lp[lower.tri(Lp, diag = TRUE)] <- par[1:(d * (d + 1) / 2)]
  Le <- matrix(0, d, d); Le[lower.tri(Le, diag = TRUE)] <- par[-(1:(d * (d + 1) / 2))]
  Psi <- Lp %*% t(Lp); E <- Le %*% t(Le)
  Sig <- make_sigma(R, Psi, E)
  ch <- tryCatch(chol(Sig), error = function(e) NULL)
  if (is.null(ch)) return(1e10)
  logdet <- 2 * sum(log(diag(ch)))
  # -2/n loglik (drop const): logdet + tr(Sig^{-1} S)
  logdet + sum(diag(chol2inv(ch) %*% S))
}
p0 <- c(t(chol(pj$Psi))[lower.tri(pj$Psi, diag = TRUE)],
        t(chol(pj$E))[lower.tri(pj$E, diag = TRUE)])
opt <- optim(p0, negll, method = "BFGS", control = list(reltol = 1e-12, maxit = 500))
Lp <- matrix(0, d, d); Lp[lower.tri(Lp, diag = TRUE)] <- opt$par[1:3]; Psi_mle <- Lp %*% t(Lp)
Le <- matrix(0, d, d); Le[lower.tri(Le, diag = TRUE)] <- opt$par[4:6]; E_mle <- Le %*% t(Le)
cat("max |Psi_blockavg - Psi_MLE| =", max(abs(pj$Psi - Psi_mle)), "\n")
cat("max |E_blockavg   - E_MLE|   =", max(abs(pj$E - E_mle)), "\n")
cat("negll(blockavg) - negll(MLE) =", negll(p0) - opt$value, " (>= 0 expected, ~0)\n\n")

# ---- (2) coefficient consistency under MCAR ----------------------------------
cat("=== (2) vector Conger/Fleiss/BP consistency (block-CS truth) ===\n")
W <- matrix(c(1, .3, .3, 1.2), 2, 2)
# per-feature continuous BP chance on a working range; D_e^BP = sum_a W_aa * range_a^2/6
ranges <- c(6, 6)
bp_chance <- sum(diag(W) * ranges^2 / 6)
truth <- vec_coef(mu_true, Sig_true, R, d, W, bp_chance)
big <- vec_coef(mu_hat, make_sigma(R, pj$Psi, pj$E), R, d, W, bp_chance)
cat("truth         :", round(truth, 4), "\n")
cat("block-CS @n=8k:", round(big, 4), "\n")
cat("max |diff|    :", round(max(abs(truth - big)), 4), "\n\n")

# d=1 reduction sanity: block-CS coincides with the scalar t1/t2/t3 forms
cat("=== d=1 reduction: vec_coef == scalar t1/t2/t3 forms ===\n")
R1 <- 5; mu1 <- c(0, .3, .6, .9, 1.2)
Psi1 <- matrix(1.0); E1 <- matrix(.5)
S1 <- make_sigma(R1, Psi1, E1)
W1 <- matrix(1)
vc <- vec_coef(mu1, S1, R1, 1, W1, bp_chance = 36 / 6)
t1 <- sum(S1); t2 <- sum(diag(S1)); t3 <- sum((mu1 - mean(mu1))^2)
sc_conger <- (t1 - t2) / ((R1 - 1) * t2 + R1 * t3)
sc_fleiss <- (t1 - t2 - t3) / ((R1 - 1) * (t2 + t3))
cat("vec  :", round(vc[1:2], 6), "\n")
cat("scalar:", round(c(sc_conger, sc_fleiss), 6), "\n")
cat("max |diff| =", max(abs(vc[1:2] - c(sc_conger, sc_fleiss))), "\n")
