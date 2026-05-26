#' Agreement study from Zapf et. al (2016)
#'
#' Agreement study (n = 200) from Zapf et al. (2016) in wide format. There are
#'   `50` items `4` judges, and ratings from `1` to `5`. It is
#'    the case that the same set four judges rated every item, hence this
#'    data is suitable for Cohen's kappa.
#'
#' @usage dat.zapf2016
#'
#' @format A `n` times `R` matrix. There are `n = 50` row corresponding to the
#'   different items. Each of the `R = 4` columns contains the ratings of
#'   the `j`th judge.
#'
#' @keywords datasets
#'
#' @references
#' Zapf, A., Castell, S., Morawietz, L. et al. Measuring inter-rater
#' reliability for nominal data <U+2013> which coefficients and confidence
#' intervals are appropriate?. BMC Med Res Methodol 16, 93 (2016).
#' https://doi.org/10.1186/s12874-016-0200-9
#'
#' @source
#' <https://biomedcentral.com/articles/10.1186/s12874-016-0200-9#Sec14>
"dat.zapf2016"
