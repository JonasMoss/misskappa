# misskappa repo-root recipes.
#
# C++ library lives at the repo root (include/, src/, tests/).
# R bindings live under r-package/.
# Manuscripts live under papers/<slug>/, each with its own justfile.

set shell := ["bash", "-euo", "pipefail", "-c"]

default:
  @just --list

# C++ build + tests (dev preset: debug, no sanitizers).
build:
  @cmake --preset dev
  @cmake --build --preset dev --parallel $(nproc)

test: build
  @ctest --preset dev

fmt-cpp:
  @formatter="$(command -v clang-format || command -v clang-format-18 || command -v clang-format-17 || true)"; \
    if [[ -z "$formatter" ]]; then \
      echo "clang-format not found; install clang-format or clang-format-18."; \
      exit 127; \
    fi; \
    rg --files include src tests/unit r-package/src \
      | rg '\.(c|cc|cpp|cxx|h|hpp)$' \
      | xargs "$formatter" -i

# Release build (no sanitizers). Used by the C++ unit tests (ctest).
opt:
  @cmake --preset opt
  @cmake --build --preset opt --parallel $(nproc)

test-opt: opt
  @ctest --preset opt

# C++ test coverage (clang source-based instrumentation + llvm-cov).
# Prints a per-file summary; full HTML lands in build-cov/coverage-html/.
cpp-cov:
  @dev/cpp-coverage.sh

# R package test coverage (covr over the testthat suite). Vendors first so the
# package's C++ is fresh; covr also reports that C++ as hit through Rcpp.
r-cov: vendor
  @Rscript dev/r-coverage.R

# Both coverage suites.
cov: cpp-cov r-cov

# Vendor the canonical C++ (src/ + include/) into r-package/src/ so the package
# is self-contained. Run after any C++ change — the R recipes below depend on it,
# so `just r-install` / `r-check` always build against fresh-vendored sources.
vendor:
  @dev/vendor-cpp.sh

# Fail if the vendored copies are stale (canonical changed without re-vendoring,
# a vendored copy was hand-edited, or a new canonical source wasn't vendored).
# For CI / pre-commit.
vendor-check: vendor
  @if [ -n "$(git status --porcelain -- r-package/src)" ]; then \
     echo "r-package/src is out of sync — run \`just vendor\` and commit:"; \
     git status --porcelain -- r-package/src; \
     exit 1; \
   fi

# Reinstall the R package against fresh-vendored sources.
r-docs:
  @cd r-package && Rscript -e 'roxygen2::roxygenise()'

r-install: vendor r-docs
  @MAKEFLAGS="-j$(nproc)" R CMD INSTALL --preclean r-package

# Fast DEV install — for day-to-day C++/logic iteration, NOT release. Compiles
# only the R glue (RcppExports.cpp + rcpp_glue.cpp) and links it against the
# opt-preset libmisskappa.a, so an install does a lib-only cmake build + a 2-file
# compile instead of recompiling the whole vendored library. Builds from a
# throwaway mirror (build-rdev/) with the dev-only dev/r-makevars-dev swapped in
# for src/Makevars, so the committed, CRAN-self-contained r-package/src/Makevars
# is never touched.
#
# Deliberately does NOT run r-docs: roxygenise() calls pkgbuild::compile_dll(),
# which recompiles all 14 files in r-package/ (the very thing we're avoiding).
# When you change the R-facing interface — add/remove a [[Rcpp::export]],
# edit roxygen/@export, touch NAMESPACE — run `just r-install` (or `just r-docs`)
# to regenerate RcppExports/NAMESPACE/man. r-dev only re-links existing exports
# against your new C++.
r-dev:
  @cmake --preset opt
  @cmake --build --preset opt --target misskappa --parallel $(nproc)   # lib only — skip the unit-test exe
  @root="$(pwd)"; \
    rsync -a --delete \
      --exclude='*.o' --exclude='*.so' --exclude='*.dll' --exclude='symbols.rds' \
      r-package/ build-rdev/; \
    cp dev/r-makevars-dev build-rdev/src/Makevars; \
    MISSKAPPA_ROOT="$root" MAKEFLAGS="-j$(nproc)" R CMD INSTALL build-rdev

r-check: vendor r-docs
  @cleanup() { rm -f misskappa_*.tar.gz; }; \
    cleanup; \
    trap cleanup EXIT; \
    R CMD build r-package; \
    R CMD check --no-manual misskappa_*.tar.gz

# Regenerate irrCAC oracle fixtures (requires R + irrCAC installed).
regen-oracle:
  @Rscript tests/tools/regen_oracle.R

# Build pkgdown site. Requires the R packages pkgdown and quarto; vendors first
# because pkgdown installs/loads the R package. If Pandoc is only available via
# Quarto, expose Quarto's bundled copy to rmarkdown/pkgdown.
docs-r: vendor
  @cd r-package && \
    if [[ -z "${RSTUDIO_PANDOC:-}" ]] && ! command -v pandoc >/dev/null 2>&1 && command -v quarto >/dev/null 2>&1; then \
      quarto_bin="$(readlink -f "$(command -v quarto)")"; \
      quarto_root="$(dirname "$(dirname "$quarto_bin")")"; \
      quarto_pandoc="$quarto_root/bin/tools/x86_64"; \
      if [[ -x "$quarto_pandoc/pandoc" ]]; then \
        export RSTUDIO_PANDOC="$quarto_pandoc"; \
      fi; \
    fi; \
    Rscript -e 'roxygen2::roxygenise(); pkgdown::build_site()'

# Build C++ API reference into the pkgdown output tree. Requires doxygen.
docs-cpp:
  @doxygen docs/Doxyfile

# Build the combined local documentation site under docs/site/.
docs: docs-r docs-cpp

docs-clean:
  @rm -rf docs/site docs/doxygen

# Render README.qmd -> README.md against the freshly built package.
readme: r-install
  @quarto render README.qmd

# Delegate to a paper-local justfile. Slug is one of: combined, ipw, fiml, quadratic.
# Example: just paper ipw pdf
paper slug *args:
  @cd papers/{{slug}} && just {{args}}

clean:
  @rm -rf build-dev build-opt build-rdev
  @rm -rf *.Rcheck
  @rm -f misskappa_*.tar.gz
  @rm -f r-package/src/*.o r-package/src/*.so r-package/src/*.dll r-package/src/symbols.rds
