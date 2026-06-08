#' Fleiss (1971) psychiatric diagnoses dataset
#'
#' Counts-format rating data from Fleiss (1971), n = 30 subjects, six raters,
#' five categories. Each row is one subject; each column counts the number of
#' raters who assigned that subject to that category.
#'
#' @format A 30 x 5 integer matrix.
#' @source Fleiss, J. L. (1971). Measuring nominal scale agreement among many
#'   raters. \emph{Psychological Bulletin}, 76(5), 378-382.
#' @examples
#' kappa_counts(dat.fleiss1971, estimator = "pairwise")
#' @keywords datasets
"dat.fleiss1971"

#' Gwet (2014) example ratings with missing data
#'
#' Raw rating data used as an example in Gwet (2014) and carried over from the
#' legacy `misskappa` package: 20 subjects, five raters, four categories
#' (`0`--`3`). Some ratings are missing (`NA`), so it is a compact worked
#' example for the missing-data categorical estimators.
#'
#' @format A 20 x 5 data frame; columns `rater1`--`rater5` hold category codes
#'   `0`--`3`, with `NA` marking unobserved ratings.
#' @source Carried over from the legacy `misskappa` package; appears as an
#'   example in Gwet, K. L. (2014), \emph{Handbook of Inter-Rater Reliability}.
#' @examples
#' kappa(dat.gwet2014, estimator = "ipw", weight = "linear")
#' @keywords datasets
"dat.gwet2014"

#' Klein (2018) example ratings with missing data
#'
#' Raw rating data carried over from the legacy `misskappa` package: 10
#' subjects, five raters, three categories (`1`--`3`), with some ratings missing
#' (`NA`). A minimal example for the missing-data categorical estimators.
#'
#' @format A 10 x 5 data frame; columns `rater1`--`rater5` hold category codes
#'   `1`--`3`, with `NA` marking unobserved ratings.
#' @source Carried over from the legacy `misskappa` package (Klein, 2018).
#' @references Klein, D. (2018). Implementing a general framework for assessing
#'   interrater agreement in Stata. \emph{The Stata Journal}, 18(4), 871-901.
#'   \doi{10.1177/1536867X1801800408}
#' @examples
#' kappa(dat.klein2018, estimator = "ipw")
#' @keywords datasets
"dat.klein2018"

#' Zapf (2016) example ratings
#'
#' Complete raw rating data carried over from the legacy `misskappa` package:
#' 50 subjects rated by four raters on a five-point scale (`1`--`5`). Having no
#' missing entries, it is a convenient base for the complete-data coefficients
#' and for illustrating the equality tests across rater pairs.
#'
#' @format A 50 x 4 data frame; columns `Rater A`--`Rater D` hold category
#'   codes `1`--`5`.
#' @source Carried over from the legacy `misskappa` package (Zapf et al., 2016).
#' @examples
#' kappa(dat.zapf2016, estimator = "pairwise")
#' @keywords datasets
"dat.zapf2016"

#' McDuff & Girard (2019) smiling-images agreement study
#'
#' Agreement study from McDuff & Girard (2019) in long format: MTurk judges
#' rate images of people smiling. The design is large and incomplete --- many
#' items and judges, and no judge rates every image --- which makes it a useful
#' missing-data example. There are 273 items, 121 judges, and 1365 rows.
#'
#' @format A 1365 x 4 data frame with one row per (item, judge) rating:
#'   \describe{
#'     \item{rating_positive}{numeric; rated positivity of the image.}
#'     \item{rating_smile}{numeric; rated presence of a smile in the image.}
#'     \item{item}{integer; the image rated.}
#'     \item{judge}{integer; the MTurk judge giving the rating.}
#'   }
#' @source The file \code{04_MTurk.csv} at \url{https://osf.io/n4grd/}.
#' @references McDuff, D., & Girard, J. M. (2019). Democratizing Psychological
#'   Insights from Analysis of Nonverbal Behavior. \emph{2019 8th International
#'   Conference on Affective Computing and Intelligent Interaction (ACII)},
#'   220-226. \doi{10.1109/ACII.2019.8925503}
#' @examples
#' # Reshape the long table to an (items x judges) matrix for one attribute,
#' # then estimate the scored agreement coefficient on the incomplete grid.
#' smile <- with(dat.mcduff2019,
#'               tapply(rating_smile, list(item, judge), function(z) z[1]))
#' smile <- matrix(as.numeric(smile), nrow = nrow(smile))
#' kappa(smile, estimator = "pairwise")
#' @keywords datasets
"dat.mcduff2019"

#' Holzinger & Swineford (1939) mental-ability battery
#'
#' The classic factor-analysis teaching dataset: 301 seventh- and eighth-grade
#' students from two schools each take nine cognitive-ability tests. The nine
#' continuous scored items split into three subscales --- *visual* (`x1`--`x3`),
#' *textual* (`x4`--`x6`), and *speed* (`x7`--`x9`) --- which makes it a compact,
#' complete real battery for coefficient alpha: estimate a subscale's alpha with
#' [alpha()], compare reliabilities across the three subscales with a G-way
#' [alpha_test()], or compare a subscale across the two schools with an
#' independent-sample [alpha_test()]. Being complete, it is also a convenient
#' base for illustrating the missing-data estimators under deliberate amputation.
#' The data are factual 1939 measurements and are in the public domain.
#'
#' @format A 301 x 15 data frame:
#'   \describe{
#'     \item{id}{subject identifier.}
#'     \item{sex}{1 = male, 2 = female.}
#'     \item{ageyr, agemo}{age, years and months.}
#'     \item{school}{factor; `"Pasteur"` (n = 156) or `"Grant-White"` (n = 145).}
#'     \item{grade}{school grade.}
#'     \item{x1, x2, x3}{visual subscale: visual perception, cubes, lozenges.}
#'     \item{x4, x5, x6}{textual subscale: paragraph comprehension, sentence
#'       completion, word meaning.}
#'     \item{x7, x8, x9}{speed subscale: speeded addition, speeded dot counting,
#'       speeded discrimination of straight and curved capitals.}
#'   }
#' @source The `HolzingerSwineford1939` data frame distributed with the
#'   \pkg{lavaan} package; originally Holzinger & Swineford (1939).
#' @references
#'   Holzinger, K. J., & Swineford, F. (1939). \emph{A study in factor analysis:
#'   The stability of a bi-factor solution}. Supplementary Educational
#'   Monographs, No. 48. University of Chicago.
#'
#'   Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
#'   \emph{Journal of Statistical Software}, 48(2), 1-36.
#'   \doi{10.18637/jss.v048.i02}
#' @examples
#' # Coefficient alpha for the textual subscale (continuous items).
#' alpha(as.matrix(dat.holzinger1939[, c("x4", "x5", "x6")]), estimator = "nt_fiml")
#' @keywords datasets
"dat.holzinger1939"

#' Vanbelle (2019) crackles lung-sound agreement study (vector-valued ratings)
#'
#' The CRACKLES auscultation study: 20 patients from the Tromso cohort, each
#' recorded at six chest sites, classified for the presence of crackles
#' (`0`/`1`) by 28 observers from seven groups of four (international experts
#' `EXP`, general practitioners from Norway `NOR`, Russia `RUS`, Wales `WAL` and
#' the Netherlands `NLD`, pulmonologists `PUL`, and medical students `STU`).
#' Because each observer rates a patient at all six sites, every rating is a
#' six-component *vector*, which makes this the worked example for the
#' vector-valued (component-separable) path of [kappa()]: pass the array and
#' pick a component loss with `weight`. Ratings are complete (no missing
#' entries). The data are factual binary classifications and are treated as
#' public domain.
#'
#' @format A 20 x 28 x 6 numeric array of `0`/`1` crackle classifications with
#'   named dimensions:
#'   \describe{
#'     \item{patient}{20 patients (`"1"`--`"20"`).}
#'     \item{rater}{28 observers, four per group: `EXP1`--`EXP4`, `NOR1`--`NOR4`,
#'       `RUS1`--`RUS4`, `WAL1`--`WAL4`, `NLD1`--`NLD4`, `PUL1`--`PUL4`,
#'       `STU1`--`STU4`.}
#'     \item{site}{six chest sites: `U1`, `U2` (upper posterior, left/right),
#'       `L1`, `L2` (lower posterior), `A1`, `A2` (anterior).}
#'   }
#' @source The `CRACKLES` data frame distributed with the \pkg{multiagree}
#'   package (Vanbelle), reshaped into a subjects-by-raters-by-features array.
#' @references
#'   Vanbelle, S. (2019). Asymptotic variability of (multilevel) multirater
#'   kappa coefficients. \emph{Statistical Methods in Medical Research},
#'   28(10-11), 3012-3026. \doi{10.1177/0962280218794733}
#'
#'   Aviles-Solis, J. C., Vanbelle, S., Halvorsen, P. A., et al. (2017).
#'   International perception of lung sounds: a comparison of classification
#'   across some European borders. \emph{BMJ Open Respiratory Research}, 4(1),
#'   e000250. \doi{10.1136/bmjresp-2017-000250}
#' @examples
#' # Vector-valued agreement across the six sites for the four expert observers,
#' # with the default Hamming component loss (count of disagreeing sites).
#' kappa(dat.vanbelle2019[, 1:4, ], estimator = "pairwise")
#' @keywords datasets
"dat.vanbelle2019"
