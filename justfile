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
  @cmake --build --preset dev

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

# Release build (no sanitizers). This is what the R package links against.
opt:
  @cmake --preset opt
  @cmake --build --preset opt

test-opt: opt
  @ctest --preset opt

# Reinstall the R package against a fresh opt build.
r-install: opt
  @R CMD INSTALL --preclean r-package

irrcacsmoke-install:
  @R CMD INSTALL --preclean dev/irrcacsmoke

r-check: opt irrcacsmoke-install
  @cleanup() { rm -f misskappa_*.tar.gz; }; \
    cleanup; \
    trap cleanup EXIT; \
    R CMD build r-package; \
    MISSKAPPA_INCLUDE="$PWD/include" \
    MISSKAPPA_LIB="$PWD/build-opt/libmisskappa.a" \
    R CMD check --no-manual misskappa_*.tar.gz

# Regenerate irrCAC oracle fixtures (requires R + irrCAC installed).
regen-oracle:
  @Rscript tests/tools/regen_oracle.R

# Build pkgdown site. Requires the R packages pkgdown and quarto; builds the
# opt library first because pkgdown installs/loads the R package.
docs-r: opt
  @cd r-package && Rscript -e 'roxygen2::roxygenise(); pkgdown::build_site()'

# Build C++ API reference into the pkgdown output tree. Requires doxygen.
docs-cpp:
  @doxygen docs/Doxyfile

# Build the combined local documentation site under docs/site/.
docs: docs-r docs-cpp

docs-clean:
  @rm -rf docs/site docs/doxygen

# Delegate to a paper-local justfile. Slug is one of: combined, ipw, fiml, quadratic.
# Example: just paper ipw pdf
paper slug *args:
  @cd papers/{{slug}} && just {{args}}

clean:
  @rm -rf build-dev build-opt
  @rm -rf *.Rcheck
  @rm -f misskappa_*.tar.gz
  @rm -f r-package/src/*.o r-package/src/*.so r-package/src/*.dll r-package/src/symbols.rds
