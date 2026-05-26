
# quadagree <img src="man/figures/logo.png" align="right" width="163" height="85" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/quadagree.png)](https://cran.r-project.org/package=quadagree)
[![R-CMD-check](https://github.com/JonasMoss/quadagree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JonasMoss/quadagree/actions/workflows/R-CMD-check.yaml)

An `R` package for calculating and doing inference the quadratically
weighted multi-rater measures of agreement. Fleiss’ kappa, Cohen’s kappa
(Conger’s kappa), and the Brennan-Prediger coefficient. Has support for
missing values using the methods of Moss and van Oest (work in progress)
and Moss (work in progress).

## Installation

The package is not available on `CRAN` yet, so use the following command
from inside `R`:

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/quadagree")
```

## Usage

Call the `library` function and load the data of Zapf et al. (2016):

``` r
library("quadagree")
head(dat.zapf2016)
#>   Rater A Rater B Rater C Rater D
#> 1       5       5       4       5
#> 2       1       1       1       1
#> 3       5       5       5       5
#> 4       1       3       3       3
#> 5       5       5       5       5
#> 6       1       1       1       1
```

Then calculate an asymptotically distribution-free confidence interval
for ![\kappa](https://latex.codecogs.com/svg.latex?%5Ckappa "\kappa"),

``` r
fleissci(dat.zapf2016)
#> Call: fleissci(x = dat.zapf2016)
#> 
#> 95% confidence interval (n = 50).
#>     0.025     0.975 
#> 0.8418042 0.9549730 
#> 
#> Sample estimates.
#>     kappa        sd 
#> 0.8983886 0.1971016
```

You can also calculate confidence intervals for Conger’s kappa (Cohen’s
kappa) and the Brennan-Prediger coefficient.

``` r
congerci(dat.zapf2016)
#> Call: congerci(x = dat.zapf2016)
#> 
#> 95% confidence interval (n = 50).
#>     0.025     0.975 
#> 0.8430854 0.9538547 
#> 
#> Sample estimates.
#>     kappa        sd 
#> 0.8984700 0.1929226
```

## Support for missing values

The inferential methods have support for missing values, using pairwise
available information in the biased sample covariance matrix. We use the
asymptotic method of van Praag (1985).

The data from Klein (2018) contains missing values.

``` r
head(dat.klein2018)
#>   rater1 rater2 rater3 rater4 rater5
#> 1      1      2      2     NA      2
#> 2      1      1      3      3      3
#> 3      3      3      3      3      3
#> 4      1      1      1      1      3
#> 5      1      1      1      3      3
#> 6      1      2      2      2      2
```

The estimates returned by `congerci`, `quadagree` and `bpci` are
consistent.

``` r
congerci(dat.klein2018)
#> Call: congerci(x = dat.klein2018)
#> 
#> 95% confidence interval (n = 10).
#>       0.025       0.975 
#> -0.03263703  0.58005580 
#> 
#> Sample estimates.
#>     kappa        sd 
#> 0.2737094 0.4062667
```

## Supported inferential techniques

`quadagree` supports three basic asymptotic confidence interval
constructions. The asymptotically distribution-free interval, the
pseudo-elliptical interval, and the normal method.

| Method       | Description                                                                                                                                                                                                                                                                                        |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `adf`        | The asymptotic distribution free method. The method is asymptotically correct, but has poor small-sample performance.                                                                                                                                                                              |
| `elliptical` | The elliptical or pseudo-elliptical kurtosis correction. Uses the unbiased sample estimator of the common kurtosis (Joanes, 1998). Has better small-sample performance than `adf` and `normal` if the kurtosis is large and ![n](https://latex.codecogs.com/svg.latex?n "n") is small.             |
| `normal`     | Assumes normality of ![X](https://latex.codecogs.com/svg.latex?X "X"). This method is not recommended since it yields too short confidence intervals when the excess kurtosis of ![X](https://latex.codecogs.com/svg.latex?X "X") is larger than ![0](https://latex.codecogs.com/svg.latex?0 "0"). |

In addition, you may transform the intervals using one of four
transforms:

1.  The [Fisher
    transform](https://en.wikipedia.org/wiki/Fisher_transformation), or
    ![\kappa\mapsto \operatorname{artanh}(\kappa)](https://latex.codecogs.com/svg.latex?%5Ckappa%5Cmapsto%20%5Coperatorname%7Bartanh%7D%28%5Ckappa%29 "\kappa\mapsto \operatorname{artanh}(\kappa)").
    Famously used in inference for the correlation coefficient.
2.  The ![\log](https://latex.codecogs.com/svg.latex?%5Clog "\log")
    transform, where
    ![\kappa \mapsto \log(1-\kappa)](https://latex.codecogs.com/svg.latex?%5Ckappa%20%5Cmapsto%20%5Clog%281-%5Ckappa%29 "\kappa \mapsto \log(1-\kappa)").
    This is an asymptotic pivot under the elliptical model with parallel
    items.
3.  The identity transform. The default option.
4.  The
    [![\arcsin](https://latex.codecogs.com/svg.latex?%5Carcsin "\arcsin")
    transform](https://en.wikipedia.org/wiki/Inverse_trigonometric_functions).
    This transform might fail when
    ![n](https://latex.codecogs.com/svg.latex?n "n") is small, as
    negative values for
    ![\hat{\kappa}](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Ckappa%7D "\hat{\kappa}")
    is possible, but
    ![\arcsin](https://latex.codecogs.com/svg.latex?%5Carcsin "\arcsin")
    do not accept them,

The option `bootstrap` does studentized bootstrapping Efron, B. (1987)
with `n_reps` repetitions. If `bootstrap = FALSE`, an ordinary normal
approximation will be used. The studentized bootstrap intervals are is a
second-order correct, so its confidence intervals will be better than
the normal approximation when
![n](https://latex.codecogs.com/svg.latex?n "n") is sufficiently large.

## Data on wide form

Some agreement data is recorded on *wide form* instead of *long form*.
Here each row contains all the possible ratings of an item along with
the total number of ratings for that item. The data of Fleiss (1971) is
on this form

``` r
head(dat.fleiss1971)
#>   depression personality disorder schizophrenia neurosis other
#> 1          0                    0             0        6     0
#> 2          0                    3             0        0     3
#> 3          0                    1             4        0     1
#> 4          0                    0             0        0     6
#> 5          0                    3             0        3     0
#> 6          2                    0             4        0     0
```

Provided the raters are exchangeable in the sense that the ratings are
conditionally independent given the item, consistent inference for the
Fleiss’ kappa and the Brennan–Prediger coefficient is possible using
`fleiss_aggr` and `bp_aggr`.

``` r
fleissci_aggr(dat.fleiss1971)
#> Call: fleissci_aggr(x = dat.fleiss1971)
#> 
#> 95% confidence interval (n = 30).
#>      0.025      0.975 
#> 0.05668483 0.51145967 
#> 
#> Sample estimates.
#> kappa.xtx        sd 
#> 0.2840722 0.5987194
```

The results agree with `irrCAC`.

``` r
irrCAC::fleiss.kappa.dist(dat.fleiss1971, weights = "quadratic")
#>      coeff.name     coeff    stderr      conf.int    p.value        pa
#> 1 Fleiss' Kappa 0.2840722 0.1111794 (0.057,0.511) 0.01612517 0.8334722
#>          pe
#> 1 0.7673958
```

## Similar software

There are several `R` packages that calculate agreement coefficients.
The most feature complete `irrCAC`, which supports calculation and
inference for agreement coefficients with more weightings than the
quadratic. However, it does not support consistent inference in the
presence of missing data, as demonstrated in the consistency vignette.

## How to Contribute or Get Help

If you encounter a bug, have a feature request or need some help, open a
[Github issue](https://github.com/JonasMoss/quadagree/issues). Create a
pull requests to contribute. This project follows a [contributer code of
conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/code_of_conduct.md).

## References

- Moss, van Oest (work in progress). Inference for quadratically
  weighted multi-rater kappas with missing raters.

- Moss (work in progress). On the Brennan–Prediger coefficients.

- Cohen, J. (1968). Weighted kappa: Nominal scale agreement with
  provision for scaled disagreement or partial credit. Psychological
  Bulletin, 70(4), 213–220. https://doi.org/10.1037/h0026256

- Conger, A. J. (1980). Integration and generalization of kappas for
  multiple raters. Psychological Bulletin, 88(2), 322–328.
  https://doi.org/10.1037/0033-2909.88.2.322

- Fleiss, J. L. (1975). Measuring agreement between two judges on the
  presence or absence of a trait. Biometrics, 31(3), 651–659.
  https://www.ncbi.nlm.nih.gov/pubmed/1174623

- Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample
  skewness and kurtosis. Journal of the Royal Statistical Society:
  Series D (The Statistician), 47(1), 183-189.
  https://doi.org/10.1111/1467-9884.00122

- Lin, L. I. (1989). A concordance correlation coefficient to evaluate
  reproducibility. Biometrics, 45(1), 255–268.
  https://www.ncbi.nlm.nih.gov/pubmed/2720055

- Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample
  skewness and kurtosis. Journal of the Royal Statistical Society:
  Series D (The Statistician), 47(1), 183-189.
  https://doi.org/10.1111/1467-9884.00122

- Klein, D. (2018). Implementing a General Framework for Assessing
  Interrater Agreement in Stata. The Stata Journal, 18(4), 871–901.
  https://doi.org/10.1177/1536867X1801800408

- Van Praag, B. M. S., Dijkstra, T. K., & Van Velzen, J. (1985).
  Least-squares theory based on general distributional assumptions with
  an application to the incomplete observations problem. Psychometrika,
  50(1), 25–36. https://doi.org/10.1007/BF02294145

- Zapf, A., Castell, S., Morawietz, L. et al. Measuring inter-rater
  reliability for nominal data – which coefficients and confidence
  intervals are appropriate?. BMC Med Res Methodol 16, 93 (2016).
  https://doi.org/10.1186/s12874-016-0200-9
