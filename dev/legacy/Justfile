# Repository workflows
#
# Requirements:
# - `just` installed: https://github.com/casey/just
# - R installed
#
# Package lives in `code/misskappa`.

set shell := ["bash", "-euo", "pipefail", "-c"]

pkg := "code/misskappa"

# List available recipes.
default:
  @just --list

# Remove common build/check artifacts.
clean:
  @rm -rf *.Rcheck
  @rm -f misskappa_*.tar.gz
  @rm -f {{pkg}}/src/*.o {{pkg}}/src/*.so {{pkg}}/src/*.dll {{pkg}}/src/symbols.rds

# Build a source tarball into the repo root.
build:
  @R CMD build {{pkg}}

# Install the package from the source tarball.
install: build
  @R CMD INSTALL --preclean misskappa_*.tar.gz

# Run the testthat tests (requires `testthat` installed).
test:
  @R -q -e 'testthat::test_local("{{pkg}}")'

# Run `R CMD check` on the built tarball (mirrors CRAN workflow more closely).
check:
  @just clean
  @just build
  @R CMD check --no-manual misskappa_*.tar.gz

# Slightly stricter check settings.
check-cran:
  @just clean
  @just build
  @R CMD check --as-cran --no-manual misskappa_*.tar.gz
