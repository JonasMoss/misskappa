### ===========================================================================
###
###   Calculating the chi square test for two data sets.
###
### ===========================================================================
source("R/functions.R")

## The table from Gwet (2008)
tab <- matrix(c(123, 2, 3, 0), 2, 2)

perreault_leigh_test(tab)
mgdm_test(tab)

## The Zapf data.
load("R/zapf.rds")

tab <- Rfast::Table(
  x = fleissci::dat.zapf2016[, 1],
  y = fleissci::dat.zapf2016[, 2],
  names = FALSE
)

perreault_leigh_test(tab)
mgdm_test(tab)
