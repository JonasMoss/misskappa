## Missing data kappas

* **Quadratic kappa** 
  * **Journal:** Statistics in Medicine.
  * Compare to Lin's concordance coefficient,
  * Multiple raters such as Zapf 2016,
  * Missing data though only MCAR (?),
  * Correct inference without normality assumption.

* **IPW and available-case kappas**
  * **Journal:** Psychological methods
  * **Pitch:** The correct way to estimate kappas with missing data.
  * When does available-case work?
  * Why is the classical standard error wrong?
  * Emphasize "intuitive" structure of the estimator.
  * Simulation studies and examples.
  * Not much math in the main text except definitions and heuristics.
  * Talk about what the assumptions (PMCAR) actually mean.
  
* **Maximum likelihood for kappas**
 * **Journal:** Psychological methods
 * More straight-forward math.
 * Emphasize conceptual simplicity.
 * Emphasize efficiency.
 * Sims can show no big gain in efficiency compared to IPW.
 * Comparison to IPW and what to recommend.
 
* **R package**
  * R journal.
  * Comparisons to existing packages etc.
  * Emphasize computatinonal efficiency and ease of use.
  
* **Application paper?**
  * Some kind of business journal collecting the results nicely.
  * Or even a low-tier ML journal?
  * Many options here.
  * Grand simulation paper with master student perhaps?
  

## Future papers?


* **Agreement models: Simuation and interpretation**
  * Discrete data
    * Skill-difficulty models
    * More general models
  * Continuous data

* **Anti-Gwet**
  * Show AC1 / AC2 "simplifications" are not needed.
  * The justification for AC1 is *wrong* and the generalization does not follow.
  * Something on "paradox of kappa".
  * **Constructive:**
    Show how to do inference estimation for the **correct** AC1/AC2 with U-stats,
    or ML. Along with the incorrect ones. (This is a minor modification to the 
    the code we already have, just different denominators. Adding it to ML
    would be trivial - but more work for U-stats. I do NOT want it in my package
    though.)

* **Anti-Krippendorff.**
  * No population estimand.
  * Superficial flexibility via arbitrary “loss functions.”
  * Naive treatment of missing data.
  * Estimator-driven, not principle-driven.
  * False impression of philosophical superiority.

* **Kappas with non-MAR missingness.**
  * Use the IPW setup with designed missingness patters.
  * Or IPW with regression.
  * Can probably find asymptotic covariance matrix.

* **Efficient estimation of kappas with continuous MAR data.**
  * "Analogue"" of ML for continuous data.
  * Loads of work, little payoff?
  * Need *something like* double/debiased machine learning.

  
  
