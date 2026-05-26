## Design doc

### **`misskappa.hpp` C++ Library: Core Design Requirements**

**Project Philosophy:** To create a high-quality, standalone, procedural C++ header-only library for calculating inter-rater agreement coefficients (Fleiss, Conger, etc.) with robust support for missing data via an EM algorithm.

**1. Portability and Dependencies:**
*   The library must be self-contained within the single `misskappa.hpp` header file.
*   **No Hard Dependencies on R/Rcpp:** The library must not call any functions from R or Rcpp. All functionality (e.g., probability distribution sampling, quantile functions) must be implemented using the C++ Standard Library.
*   **Minimal External Dependencies:** The only required dependency is a BLAS/LAPACK-enabled linear algebra library. The code is written for the **Armadillo** API. It must be compatible with a standard Armadillo installation (for portability to Python, etc.), not just `RcppArmadillo`.
    *   *Note: This implies that all types should be standard C++ or `arma::*` types. `Rcpp::*` types are forbidden within the `.hpp` file.*

**2. Programming Paradigm and Style:**
*   **Procedural, Not Object-Oriented:** The library will follow a procedural C-style approach. Functionality is organized into namespaces, not classes. Data is passed between functions using simple `structs` without methods or private members.
*   **Minimal Use of Advanced C++ Features:** Templates will be used sparingly and only where they significantly reduce code duplication without obfuscating logic (e.g., for error handling `Result<T>` and data-type generic functions).
*   **Style Guide Adherence:** The code will follow a consistent style similar to the Google C++ Style Guide: `snake_case` for functions/variables, `PascalCase` for types, and `kConstant` for enumerators and compile-time constants.

**3. Error Handling:**
*   **No C++ Exceptions:** The library must not use C++ exceptions for control flow or to signal recoverable errors.
*   **Result-Based Error Propagation:** All functions that can fail must return a `misskappa::Result<T>` struct. This struct explicitly packages either a success value (`T`) or an error status and message, forcing the caller to handle potential failures.

**4. Performance:**
*   **Efficiency is Key:** The implementation must be computationally efficient. This includes:
    *   Preprocessing data to operate on compact representations (unique pattern counts).
    *   Using memoization for expensive repeated calculations (e.g., combinations).
    *   Implementing efficient resampling strategies (e.g., resampling from a subject pool instead of the full raw data in the non-parametric bootstrap).

**5. API Stability and Extensibility:**
*   **Stable Public API:** The primary functions exposed in the `misskappa` namespace (`preprocess_*`, `run_em`, `calculate_estimate`, etc.) form the stable public API. This interface should be logical and sufficient for building external wrappers (like R or Python bindings) without modification.
*   **Extensible by Design:** The architecture must make it easy to add new agreement coefficients. The "estimator factory" pattern, which uses `std::function`, is the core of this design. Adding a new coefficient should only require writing a small, self-contained factory function.

**6. Compatibility with math.**
*   The notation and names must be compatible with our manuscript.
*   In particular, we shall only use *disagreement metrics*.

**7. RCPP / R / testing.**
*   There is a companion `misskappa.cpp` that binds to `R` and should work without additional wrappers.
*   dat.zapf2016 and dat.fleiss1971 shall be used for testing complete data against irrCAC. 

**Future Work:**

1. **Fix interface.** The package MUST NOT export options functions. These should be passed as `...` OR as incomplete options arguments that are documented extensively in the docs.
  This is important as it clobbers the interface.
2. Connect parametric bootstrap to `R`.
3. **Priority: 2/10** The hashing uses strings. It could use ranks or vector hashing instead. 


