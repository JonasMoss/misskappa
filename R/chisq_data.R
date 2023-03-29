### ===========================================================================
###
### We show that Gwet's (2008) table 4 is compatible with the Perreault--Leigh
###    model. (Is this at all informative? I don't really think so.)
###
### ===========================================================================
source("R/functions.R")

tab <- matrix(c(123, 2, 3, 0), 2, 2)
perreault_leigh_test(tab)
mgdm_test_2(tab)

tab <- Rfast::Table(
  x = fleissci::dat.zapf2016[, 1],
  y = fleissci::dat.zapf2016[, 2],
  names = FALSE
)
perreault_leigh_test(tab)
mgdm_test_2(tab)
