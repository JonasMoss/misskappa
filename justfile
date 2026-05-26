# misskappa repo-root recipes.
#
# C++ library lives at the repo root (include/, src/, tests/).
# R bindings live under r-package/.
# Manuscript lives under paper/ with its own justfile.

set shell := ["bash", "-euo", "pipefail", "-c"]

default:
  @just --list

# C++ build + tests (dev preset: asan + ubsan).
build:
  @cmake --preset dev
  @cmake --build --preset dev

test: build
  @ctest --preset dev

# Release build (no sanitizers). This is what the R package links against.
opt:
  @cmake --preset opt
  @cmake --build --preset opt

test-opt: opt
  @ctest --preset opt

# Reinstall the R package against a fresh opt build.
r-install: opt
  @R CMD INSTALL --preclean r-package

r-check: opt
  @R CMD build r-package
  @R CMD check --no-manual misskappa_*.tar.gz
  @rm -f misskappa_*.tar.gz

# Regenerate irrCAC oracle fixtures (requires R + irrCAC installed).
regen-oracle:
  @Rscript tests/tools/regen_oracle.R

# Delegate to the paper-local justfile for the manuscript.
paper *args:
  @cd paper && just {{args}}

clean:
  @rm -rf build-dev build-opt
  @rm -rf *.Rcheck
  @rm -f misskappa_*.tar.gz
  @rm -f r-package/src/*.o r-package/src/*.so r-package/src/*.dll r-package/src/symbols.rds
