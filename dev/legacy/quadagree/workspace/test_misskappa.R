weights <- c("quadratic", "ordinal", "linear", "radical", "ratio", "circular", "bipolar", "identity")

set.seed(313)
x <- as.matrix(dat.zapf2016)
#x[sample(length(x), 40)] <- NA
storage.mode(x) <- "integer"
x
em_opts <- list(tol = 1e-8, max_iter = 1000, prune_tol = 1e-12, start_alpha = 0.5)
x_fleiss <- as.matrix(dat.fleiss1971)
storage.mode(x) <- "integer"


sapply(weights, \(w) {
  fit <- conger_cpp(x, 5, weight_type = w, values = c(1,2,3,4,5), em_options = em_opts)
  fit2 <- conger(x, weights = w)
  list(weight = w,
       kappa = fit$estimate,
       kappa_irr = fit2[1],
       stderr = sqrt(fit$variance),
       stderr_irr = fit2[2])
})

sapply(weights, \(w) {
  fit <- fleiss_cpp(x_fleiss, 6, weight_type = w, values = c(1,2,3,4,5), em_options = em_opts)
  fit2 <- fleiss(x_fleiss, weights = w)
  list(weight = w,
       kappa = fit$estimate,
       kappa_irr = fit2[1],
       stderr = sqrt(fit$variance),
       stderr_irr = fit2[2])
})


sapply(weights, \(w) {
  fit <- bp_cpp(x_fleiss, 6, weight_type = w, values = c(1,2,3,4,5), em_options = em_opts, NULL)
  fit2 <- bp(x_fleiss, weights = w)
  list(weight = w,
       kappa = fit$estimate,
       kappa_irr = fit2[1],
       stderr = sqrt(fit$variance),
       stderr_irr = fit2[2])
})

sapply(weights, \(w) {
  fit <- fleiss_raw_cpp(x, 5, weight_type = w, values = c(1,2,3,4,5), em_options = em_opts, NULL)
  fit2 <- fleiss(to_counts_matrix(x), weights = w)
  list(weight = w,
       kappa = fit$estimate,
       kappa_irr = fit2[1],
       stderr = sqrt(fit$variance),
       stderr_irr = fit2[2])
})


sapply(weights, \(w) {
  fit <- bp_raw_cpp(x, 5, weight_type = w, values = c(1,2,3,4,5), em_options = em_opts, NULL)
  fit2 <- bp(to_counts_matrix(x), weights = w)
  list(weight = w,
       kappa = fit$estimate,
       kappa_irr = fit2[1],
       stderr = sqrt(fit$variance),
       stderr_irr = fit2[2])
})


analysis_opts <- list(
  conf_level = 0.95,
  ci_transform = "fisher", # "identity", "fisher", "log", "arcsin"
  bootstrap_method = "nonparametric-t", # "none", "nonparametric", "nonparametric-t"
  bootstrap_reps = 1000,
  seed = 123
)

bp_cpp(x_fleiss, 6, weight_type = "quadratic", values = c(1,2,3,4,5), em_options = em_opts, analysis_options = analysis_opts)

fleiss_cpp(counts_matrix, 4, "identity", values = c(1,2,3,4), em_options = em_opts, analysis_options = NULL)
