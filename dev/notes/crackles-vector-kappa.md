# CRACKLES vector-kappa notes

Working note for the vector-valued agreement pilot in
`experiments/26-crackles-vector-kappa/`.

Status: the active implementation target is now the component-separable vector
kappa in `dev/notes/component-separable-vector-kappa.md`. The original
CRACKLES pilot used a pure-R squared-loss mean/covariance contraction with
normal-FIML imputation as a feasibility check. The experiment should now call
the internal `kappa_vector()` wrapper for pairwise and IPW component-missing
estimators; categorical full-profile FIML is deferred.

## What the data are

`dev/dat/CRACKLES.rda` is the CRACKLES / Tromso lung-sound agreement
example used by Vanbelle (2019). It comes from the lung-sound perception
study of Aviles-Solis et al. (2017).

The unit structure is:

- 20 subjects / patients.
- 6 recordings per subject.
- The 6 recordings are three thorax locations on each side of the chest:
  upper posterior, lower posterior, and anterior.
- 28 observers, arranged as seven groups of four observers.
- Binary response: whether crackles are present in the recording.

The data frame has 120 rows because it is stored at the recording level:
`20 patients x 6 recordings = 120 rows`.

## Site coding

The site indicators encode thorax location, not missingness:

- `UP == 1`, `LO == 0`: upper posterior thorax.
- `UP == 0`, `LO == 1`: lower posterior thorax.
- `UP == 0`, `LO == 0`: anterior thorax.

Within each patient, the observed row order is:

1. upper posterior
2. upper posterior
3. lower posterior
4. lower posterior
5. anterior
6. anterior

The paired rows are the two sides of the chest. The source papers describe
the design as six body sites, or equivalently three thorax locations on each
side.

Implication: the vector pilot should describe the patient-level response as
a 6-site binary vector, but the clinically meaningful grouping is
`location x side`, with locations `U`, `L`, and `A`. The current experiment
uses the six sites as an unordered vector with `W = I`; that is fine for a
first Hamming-loss pilot, but any location-specific analysis should collapse
or stratify by `U/L/A` as Vanbelle does.

## Observer classes

The 28 rating columns are seven observer classes with four raters each:

| Prefix | Vanbelle label | Observer class |
| --- | --- | --- |
| `EXP` | `EXP` | International lung-sound experts / researchers |
| `NOR` | `NOR` | General practitioners from Norway |
| `RUS` | `RUS` | General practitioners from Russia |
| `WAL` | `WAL` | General practitioners from Wales |
| `NLD` | `NLD` | General practitioners from The Netherlands |
| `PUL` | `PLN` | Pulmonologists, University Hospital of North Norway |
| `STU` | `STU` | Sixth-year medical students, Faculty of Health Sciences in Tromso |

The `PUL` / `PLN` mismatch is only a naming mismatch between the data object
and Vanbelle's tables. It is the same pulmonologist group.

## Exchangeability and estimands

There are two natural readings of the observer columns.

The defensible primary reading is **within-class agreement**:

- Analyze `EXP1:EXP4`, `NOR1:NOR4`, etc. as seven separate four-rater panels.
- Within each class, raters are plausibly exchangeable for a multirater
  agreement coefficient.
- This matches Vanbelle's Table 3, which reports agreement by observer class
  and thorax location, plus an overall row.

The pooled 28-rater reading is only descriptive:

- Pooling all 28 observers treats experts, GPs from different countries,
  pulmonologists, and students as one rater population.
- That can be useful as a stress test or a broad "all observers" summary.
- It should not be presented as the main interchangeable-rater estimand unless
  the target population is explicitly "the heterogeneous panel of all observer
  classes."

## What the vector coefficient means here

For each patient, a rater gives a binary vector of length 6:

```text
(upper-left/right, lower-left/right, anterior-left/right)
```

With `W = I`, the vector disagreement between two raters is squared
Euclidean distance. For binary site indicators this is Hamming distance:
the number of sites where the two raters disagree about crackles.

The proposed coefficient contracts a single empirical mean/covariance pair
`(mu, Sigma)`:

- `T`: within-rater site variance, summed over raters and sites.
- `B`: same-site cross-rater covariance, summed over all rater pairs.
- `G`: between-rater spread in site-level positive-rate vectors.

This is a patient-level, vector-valued analogue of the scalar quadratic /
Frechet agreement story. It asks whether raters agree on a six-site crackle
profile, not merely whether they agree after pooling all recordings.

Important limitation: with `W = I`, cross-site covariances do not enter the
estimand directly. They can still help a missing-data estimator reconstruct
`(mu, Sigma)` if some site ratings are missing.

## Missingness status

The CRACKLES data object is complete: all 28 observers classify all 120
recordings. Missing-site behavior in the experiment is synthetic. That is
acceptable for a pilot, but it should be described as a stress test, not as
an analysis of naturally incomplete CRACKLES data.

For a paper-facing experiment:

- Keep the complete-data CRACKLES analysis as the applied illustration.
- Use synthetic deletion only to demonstrate how the moment plug-in behaves
  when sites are missing.
- If we want MAR rather than MCAR, make the deletion depend on observed
  site/location labels or observed ratings, not on hidden truth.

## References in `dev/refs`

- `Aviles-Solis et al. 2017 - International perception of lung sounds.pdf`
- `Aviles-Solis et al. 2019 - Prevalence and clinical associations of wheezes and crackles in the general population.pdf`
- `vanbelle2018.pdf`: Vanbelle (2019), direct CRACKLES / Tromso multilevel
  multirater kappa example.
- `vanbelle2017.pdf`: Vanbelle (2017), dependent kappa coefficients on
  multilevel data; relevant for comparison tests and the `multiagree` package.
