# Project description

This is a paper about agreement coefficients with missing values.

## Details

We will discuss statistical properties of two estimaters for Cohen's and Fleiss' kappa and the Brennan--Prediger coefficient. The IPW estimator and the available case estimator. Our goals are:

* to show that IPW is more efficient than available case, 
* sufficient conditions for consistency of both,
* joint inferential theory,
* cover data on two forms: where the raters are known and where the raters are unknown.
* show that Gwet's earlier inferential method is inconsistent (I do not recall the details),
* small example and small simulation.
* make an R package `misskappa` with Rcpp that calculates the stuff, 
* make a theoretical case that IPW is efficient, but not in the usual semiparametric sense. I am thinking about something like efficiency when we do not assume a well-specified missingness mechanism or something similar.

All of this should be possible to put into one paper without having it crammed provided we work well.

## Introduction and conclusion

* Eficiency through the EM algorithm and efficient missing data stuff; that will be another paper (ideally submitted just after this).
* A short literature review on missing values for kappas. There are not many papers here. 
* Mention the key papers in missing data but do not dwell on it at all. Three sentences is enough. Should cite Tsiasis, van der Vaart, and 

## Venue and size

* Psychometrika, aiming for 12 pages or so not counting appendix. 
* Online appendix is fine.
* I do not want a massive paper; we want clarity and readability first.