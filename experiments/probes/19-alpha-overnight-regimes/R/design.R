# The design: five claim-driven sub-studies (each a 2-D slice with the rest at
# defaults) plus an individually-reported torture gallery. Cells are rows; the
# overnight cost is driven by reps, not cell count.
#
#   S1  pairwise SE vs the average-n Feldt strawman          (the contribution)
#   S2  MCAR -> MAR -> MNAR consistency split
#   S3  nonnormality (IG/Pearson) & the SE/point-bias story
#   S4  off-congeneric x small-n bias
#   S5  efficiency (pairwise vs FIML) & cat-FIML cost/feasibility
#   torture  zero-overlap, t3, alpha~.97, contamination, wild item, C^p wall

CONT_EST <- "pairwise,fiml_normal,fiml_sandwich,feldt"
DISC_EST <- "pairwise,fiml_normal,fiml_sandwich,cat_fiml,feldt"

.cell <- function(substudy, label, model, p, n, alpha_target, dist, mechanism,
                  n_cat = NA_integer_, estimators = NULL, bootstrap = FALSE) {
  if (is.null(estimators)) estimators <- if (dist == "discrete") DISC_EST else CONT_EST
  data.frame(substudy = substudy, label = label, model = model, p = p, n = n,
             alpha_target = alpha_target, dist = dist, n_cat = n_cat,
             mechanism = mechanism, estimators = estimators, bootstrap = bootstrap,
             stringsAsFactors = FALSE)
}

build_design <- function(which = "all") {
  d <- list()

  # S1 -- the contribution: coverage vs n x missingness, random + planned MCAR.
  for (n in c(50, 100, 200, 500, 1000)) for (mech in c("mcar30", "mcar_planned")) {
    d[[length(d) + 1]] <- .cell("s1", sprintf("n=%d/%s", n, mech),
      "congeneric_mild", 8, n, 0.80, "normal", mech, bootstrap = (n == 200))
  }

  # S2 -- consistency split across mechanisms.
  for (mech in c("complete", "mcar30", "mar", "mnar")) {
    d[[length(d) + 1]] <- .cell("s2", mech, "congeneric_mild", 8, 300, 0.80, "normal", mech,
      bootstrap = (mech == "mcar30"))
  }

  # S3 -- nonnormality x mechanism (continuous), plus a discrete MAR companion
  # where cat-FIML is meant to rescue the normal-FIML point bias.
  for (dist in c("normal", "ig_skew", "ig_heavy", "t3", "contam")) for (mech in c("mcar30", "mar")) {
    d[[length(d) + 1]] <- .cell("s3", sprintf("%s/%s", dist, mech),
      "congeneric_mild", 8, 300, 0.80, dist, mech, bootstrap = (dist == "normal" && mech == "mcar30"))
  }
  for (mech in c("mcar30", "mar")) {
    d[[length(d) + 1]] <- .cell("s3", sprintf("discrete5/%s", mech),
      "congeneric_mild", 6, 300, 0.80, "discrete", mech, n_cat = 5L)
  }

  # S4 -- off-congeneric x small-n bias.
  for (model in c("parallel", "tau", "congeneric_mild", "congeneric_wild"))
    for (n in c(30, 50, 100, 200, 500, 1000)) {
      d[[length(d) + 1]] <- .cell("s4", sprintf("%s/n=%d", model, n),
        model, 8, n, 0.80, "normal", "mcar10")
    }

  # S5 -- efficiency sweep + cat-FIML cost/feasibility (small-p corner).
  for (mech in c("mcar10", "mcar30", "mcar50", "mcar_planned")) {
    d[[length(d) + 1]] <- .cell("s5", sprintf("ARE/%s", mech),
      "congeneric_mild", 8, 300, 0.80, "normal", mech)
  }
  for (p in c(4, 6, 8)) for (nc in c(3L, 5L)) {
    d[[length(d) + 1]] <- .cell("s5", sprintf("cost/p=%d/C=%d", p, nc),
      "congeneric_mild", p, 300, 0.80, "discrete", "mcar30", n_cat = nc)
  }

  # Torture gallery -- reported individually.
  d[[length(d) + 1]] <- .cell("torture", "zero_overlap", "congeneric_mild", 8, 400, 0.80, "normal", "zero_overlap")
  d[[length(d) + 1]] <- .cell("torture", "t3_no_4th_moment", "congeneric_mild", 8, 300, 0.80, "t3", "mcar30", bootstrap = TRUE)
  d[[length(d) + 1]] <- .cell("torture", "alpha_0.97_irregular", "parallel", 8, 200, 0.97, "normal", "mcar30")
  d[[length(d) + 1]] <- .cell("torture", "contamination", "congeneric_mild", 8, 300, 0.80, "contam", "mcar30", bootstrap = TRUE)
  d[[length(d) + 1]] <- .cell("torture", "wild_items", "congeneric_wild", 8, 200, NA_real_, "normal", "mcar30")
  d[[length(d) + 1]] <- .cell("torture", "catfiml_Cp_wall", "congeneric_mild", 12, 300, 0.80, "discrete", "mcar30",
                              n_cat = 5L, estimators = "cat_fiml")

  out <- do.call(rbind, d)
  out$cell_id <- seq_len(nrow(out))
  if (!identical(which, "all")) out <- out[out$substudy %in% which, , drop = FALSE]
  out
}
