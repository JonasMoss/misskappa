
rm(dat.fleiss1971)
dat.fleiss1971 <- as.matrix(dat.fleiss1971)
# Test with missing data
incomplete_data <- dat.fleiss1971
incomplete_data[1, ] <- c(0, 0, 0, 5, 0)
incomplete_data[2, ] <- c(0, 2, 0, 0, 2)
R_total <- 6

x <- incomplete_data
#x <- dat.fleiss1971


em_counts <- function(O, R, eps = 1e-8, maxit = 1000, sparse = TRUE) {
  state_space <- function(R, C) {
    grd <- expand.grid(rep(list(0:R), C))
    as.matrix(grd[rowSums(grd) == R, , drop = FALSE])  #  S × C
  }
  if (any(is.na(O))) O[is.na(O)] <- 0               # NA → 0

  n <- nrow(O);
  C <- ncol(O)
  S <- state_space(R, C);
  Sdim <- nrow(S)

  theta <- rep(1 / Sdim, Sdim)

  compat <- lapply(
    1:n,
    function(i) which(apply(S >= matrix(O[i, ], nrow = Sdim, ncol = C,
                                        byrow = TRUE), 1, all))
  )

  for (it in 1:maxit) {
    tau <- matrix(0, n, Sdim)
    for (i in 1:n) {
      idx        <- compat[[i]]
      w          <- theta[idx]
      tau[i,idx] <- w / sum(w)                 # E-step
    }
    new_theta <- colMeans(tau)                 # M-step
    if (max(abs(new_theta - theta)) < eps) break
    theta <- new_theta
  }

  # ---------- 3.  Louis observed information ---------------------------
  N <- colSums(tau)                            # length Sdim
  D <- Matrix::Diagonal(x = N / theta^2)

  # T_ij = Σ_i τ_{ix} τ_{iy}
  #T <- if (sparse) {
  #  Matrix::crossprod(Matrix::Matrix(tau, sparse = TRUE))
  #} else crossprod(tau)

  #ThetaOuter <- theta %*% t(theta)             # θ θᵀ
  #Info <- Matrix::crossprod(Matrix::Matrix(tau, sparse = TRUE)) /  theta %*% t(theta)


  # Louis matrix

  eps  <- .Machine$double.eps
  keep <- which(theta > eps)
  theta_pos <- theta[keep]
  tau_pos   <- tau[, keep]

  N_pos <- colSums(tau_pos)
  T     <- Matrix::crossprod(Matrix::Matrix(tau_pos, sparse = TRUE))
  Info  <- T / tcrossprod(theta_pos)          # Hessian without NaNs


  list(theta = theta,
       Info   = Info,
       tau    = tau,
       iter   = it)
}

# ---------- 4.  Example on your data -----------------------------------
R <- 6                     # 6 possible ratings per item
O <- x
eps <- 1e-8
maxit <- 1000
sparse <- TRUE


fit <- em_counts(x, R)     # x = 30×5 matrix you posted
theta_hat <- fit$theta
Info_hat  <- fit$Info

# quick sanity: rows should sum to 0 (simplex singularity)
max_row_sum <- max(abs(Matrix::rowSums(Info_hat)))
cat("max |row-sum| =", max_row_sum, "\n")   # ~1e-12 ⇒ OK



fit <- em_counts(x, R, sparse = TRUE)
theta <- fit$theta
tau   <- fit$tau

eps  <- 1e-16
keep <- which(theta > eps)
theta_pos <- theta[keep]
tau_pos   <- tau[, keep]

N_pos <- colSums(tau_pos)
T     <- Matrix::crossprod(Matrix::Matrix(tau_pos, sparse = TRUE))
Info  <- T / tcrossprod(theta_pos)          # Hessian without NaNs


## theta_pos, Info   already built without the zero cells
Spos <- length(theta_pos)
grad_f <- 2 * theta_pos           # analytic gradient
#
# # ---- A.  Trim one coordinate --------------------------
# Info_trim <- Info[-1, -1]
# cov_trim  <- solve(Info_trim)
#
# grad_trim <- grad_f[-1]           # drop same coordinate
# var_f_A   <- as.numeric( t(grad_trim) %*% cov_trim %*% grad_trim )
# SE_A      <- sqrt(var_f_A)

# ---- B.  Orthogonal projection ------------------------
library(Matrix)
One <- Matrix(1, Spos, 1)
P   <- Diagonal(Spos) - (1/Spos) * (One %*% t(One))   # centring matrix

cov_MP  <- MASS::ginv(as.matrix(Info))
cov_proj<- P %*% cov_MP %*% P                  # annihilate null space
var_f_B <- as.numeric( t(grad_f) %*% cov_proj %*% grad_f )
SE_B    <- sqrt(max(var_f_B, 0))               # guard tiny negatives
SE_B

