# Build data/dat.vanbelle2019.rda from the CRACKLES lung-sound agreement data.
# Source: the CRACKLES data frame distributed with Sophie Vanbelle's 'multiagree'
# package (CRAN); the statistical source is Vanbelle (2019) and the study design
# is Aviles-Solis et al. (2017). Saved here as vanbelle2019.csv for a
# network-free, reproducible build. The ratings are factual binary
# crackle classifications and are treated as public domain.
#
# The CSV is the raw 120 x 31 table: patient, UP, LO, and 28 observer columns.
# Each patient has six recordings in a fixed row order -- two upper-posterior,
# two lower-posterior, two anterior sites -- which we reshape into a
# subjects-by-raters-by-features array (20 patients x 28 observers x 6 sites)
# ready for kappa()'s vector-valued path. Run from the package root.

d <- read.csv("data_raw/vanbelle2019.csv", check.names = FALSE,
              stringsAsFactors = FALSE)

raters <- setdiff(names(d), c("patient", "UP", "LO"))
sites <- c("U1", "U2", "L1", "L2", "A1", "A2")

# Cross-check the fixed within-patient site order against the UP/LO indicators
# (rows: UP, UP, LO, LO, anterior, anterior).
expected_up <- c(1, 1, 0, 0, 0, 0)
expected_lo <- c(0, 0, 1, 1, 0, 0)
for (p in unique(d$patient)) {
  rows <- d[d$patient == p, ]
  if (nrow(rows) != 6L ||
      !identical(as.numeric(rows$UP), expected_up) ||
      !identical(as.numeric(rows$LO), expected_lo)) {
    stop("Unexpected site layout for patient ", p, ".")
  }
}

patients <- unique(d$patient)
dat.vanbelle2019 <- array(
  NA_real_,
  dim = c(length(patients), length(raters), length(sites)),
  dimnames = list(patient = as.character(patients),
                  rater = raters, site = sites)
)
for (i in seq_along(patients)) {
  rows <- d[d$patient == patients[i], raters, drop = FALSE]
  # rows is sites x raters; transpose to raters x sites.
  dat.vanbelle2019[i, , ] <- t(as.matrix(rows))
}
storage.mode(dat.vanbelle2019) <- "double"

save(dat.vanbelle2019, file = "data/dat.vanbelle2019.rda", compress = "xz")
rm(dat.vanbelle2019)
