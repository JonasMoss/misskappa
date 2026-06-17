# Build data/dat.tenhove2025.rda: a SYNTHETIC planned-incomplete inter-rater
# dataset for illustrating the compound-symmetry (variance-component) estimators.
#
# This is NOT real data. It is generated from a fixed-rater variance-component
# model in the *style* of the planned-incomplete ICC designs studied by
# ten Hove, Jorgensen & van der Ark (2025, "How to Estimate ICCs for Interrater
# Reliability from Planned Incomplete Data", Multivariate Behavioral Research
# 60(5)). The point of the design is that the rater co-observation graph is
# CONNECTED but NOT COMPLETE: some rater pairs never rate the same subject. That
# makes the saturated coefficients (estimator = "pairwise" / "nt_fiml" /
# "cat_fiml") non-identified, while the compound-symmetry estimators
# ("nt_fiml_cs", "pairwise_cs"), which need only a connected graph, still fit.
#
# 60 subjects, 6 raters, ordinal 1-5 ratings, three overlapping rater panels of
# three raters each:
#   subjects  1-20 : raters 1,2,3
#   subjects 21-40 : raters 3,4,5
#   subjects 41-60 : raters 5,6,1
# Co-observed pairs: {12,13,23, 34,35,45, 56,15,16}. The graph is connected
# (1-2-3-4-5-6 all linked) but pairs such as 2-4 are never co-observed.
# Run from the package root.

set.seed(2025)

n_subjects <- 60L
panels <- list(1:20, 21:40, 41:60)
panel_raters <- list(c(1L, 2L, 3L), c(3L, 4L, 5L), c(5L, 6L, 1L))

# Fixed-rater variance-component truth: y_ir = mu_r + b_i + e_ir.
mu <- c(2.9, 3.1, 3.0, 3.2, 2.8, 3.0)   # small per-rater bias (free means)
sd_b <- 1.00                            # subject (latent agreement) SD  -> psi
sd_e <- 0.55                            # residual disagreement SD       -> e

b <- stats::rnorm(n_subjects, 0, sd_b)

dat.tenhove2025 <- matrix(
  NA_real_, n_subjects, 6L,
  dimnames = list(NULL, paste0("rater", 1:6))
)
for (p in seq_along(panels)) {
  for (i in panels[[p]]) {
    for (r in panel_raters[[p]]) {
      y <- mu[r] + b[i] + stats::rnorm(1L, 0, sd_e)
      dat.tenhove2025[i, r] <- min(max(round(y), 1L), 5L)
    }
  }
}
storage.mode(dat.tenhove2025) <- "integer"

save(dat.tenhove2025, file = "data/dat.tenhove2025.rda", compress = "xz")
rm(dat.tenhove2025)
