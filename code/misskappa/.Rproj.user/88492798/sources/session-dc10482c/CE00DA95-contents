# **`misskappa` API Design Document**

## TODO.

### Extensions

* [ ] Have all functions return d, d_C, d_F, and their covariance matrix.
  * [ ] ML
  * [ ] np
  * [ ] moment
* [ ] Verify that definition of weighted estimator is incorrect! :(
* [ ] Verify the current definition of the SE for IPW is incorrect! :(

### Package internals stuff
* [x] Implement `np` for count data (might be done in the cheap way using to_counts_matrix),
  * [x] Update `test-complete-counts-equal.R`
* [x] Implement `np` and `ipw` for *continuous data*.


* [ ] Loss names and exposition?
  * [ ] Go with the mathematically defensible names?

* [ ] Method renaming:
  * "available" -> "available"
  * "quadratic" -> "moment"

* [ ] Custom weight matrices.
  * Will it play well with the weight matrix classes?
  
* [ ] Rename "kappaqp"" => "kappamom".

* [x] rename "kappa" -> "kappa_raw"

* [ ] Make simulation functionality.

### Verification stuff

* simulation.
  * **Simple.** A JSM.
  * **Complex.** A skilld-difficulty model.
* PMCAR mechanism.
  * `pmcar <- \(n, r, seed = 1)`
  * generates `pmcar` mechanism for `r` raters that isn't `MCAR`. Returns a probability tensor.

#### **Primary Estimation Functions**

These are the main functions users will call to perform the agreement analysis. They are designed to be simple and focus on the core task of estimation.

`kappa(x, method, weight, values, ...)`
`kappa_continuous(x, method, weight)`
`kappa_counts(x, method, weight, values, ...)`

**Return Value:**

All estimation functions will return a rich S3 object of class `"misskappa"`. This object is a list containing the core results necessary for all subsequent inference.

A `"misskappa"` object will contain at least:

*   `estimates`: A named numeric vector of the calculated kappa coefficients.
*   `vcov`: The full variance-covariance matrix of the estimates.
*   `ci_default`: A data frame containing the default confidence intervals (e.g., 95% Fisher's Z-transform).
*   `n_eff`: The number of subjects used in the analysis.
*   `call`: The original function call.
*   `metadata`: A list containing information about the analysis, such as the method and data type, needed by other functions.

---

#### **Core S3 Methods for the `"misskappa"` Object**

These methods provide the primary interface for users to interact with the results.

##### `print()` Method

**Purpose:** To provide a clean, concise, and immediate summary of the most important results. This is what the user sees when they type the name of the returned object.

**Usage:**
```R
fit <- kappa_raw(data, ...)
fit 
```

**Output:**
```
misskappa Agreement Analysis
Method: Maximum Likelihood 
Call: kappa(x = data, method = "ml")

Estimates and 95% Confidence Intervals (Fisher's Z-transform):
                   Estimate  Lower CI  Upper CI
Conger             0.751     0.681     0.821
Fleiss             0.749     0.679     0.819
Brennan-Prediger   0.720     0.643     0.797
```

##### `summary()` Method

**Purpose:** To perform and display detailed hypothesis tests. This is the main entry point for statistical inference.

**Usage:**
```R
s_fit <- summary(fit)
s_fit
```

**Output:**
The `summary()` method returns an object of class `"summary.misskappa"`, which has its own print method to display:
```
--- Coefficient Significance ---
(H0: Kappa = 0)
                 Estimate Std. Error z value Pr(>|z|)    
Conger             0.751      0.035      21.4   <2e-16 ***
Fleiss             0.749      0.036      20.8   <2e-16 ***
Brennan-Prediger   0.720      0.038      18.9   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

--- Comparison of Coefficients ---
(H0: Conger's Kappa - Fleiss' Kappa = 0)
Difference: 0.002, SE of Difference: 0.001, z = 2.0, p-value = 0.045
```

##### `confint()` Method

**Purpose:** To calculate and display confidence intervals with user-specified options, without re-running the core estimation.

**Usage:**
```R
# Calculate CIs using the identity transform (no transform)
confint(fit, transform = "identity")

# Calculate 99% CIs
confint(fit, level = 0.99)

# Future: Calculate bootstrap CIs (on-demand)
# confint(fit, method = "bootstrap", n.reps = 2000, type = "bca")
```

**Arguments:**
*   `object`: A `"misskappa"` object.
*   `level`: The confidence level (e.g., 0.95).
*   `transform`: The variance-stabilizing transformation (`"fisher"`, `"arcsin"`, `"identity"`).
*   `method`: The CI calculation method (`"wald"`, `"bootstrap"`).
*   `...`: Additional arguments passed to the bootstrap routine (e.g., `n.reps`, `type`).

##### `coef()` and `vcov()` Methods

**Purpose:** To provide standard, programmatic access to the core numerical results, consistent with other R model objects.

**Usage:**
```R
# Extract the named vector of estimates
coef(fit)
# >       Conger         Fleiss Brennan-Prediger 
# >    0.7511111      0.7492222        0.7203333 

# Extract the full variance-covariance matrix
vcov(fit)
# >                  Conger        Fleiss  Brennan-Prediger
# > Conger           0.001225      0.001220       0.001211
# > Fleiss           0.001220      0.001296       0.001280
# > Brennan-Prediger 0.001211      0.001280       0.001444
```

---
This API design provides a clean, layered approach. The user gets immediate, simple results from `print()`, can dig deeper with `summary()`, and has full flexibility for custom CIs and programmatic access with `confint()`, `coef()`, and `vcov()`. It's robust, idiomatic, and perfectly sets the stage for implementation.
