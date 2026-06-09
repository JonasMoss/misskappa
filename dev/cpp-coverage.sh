#!/usr/bin/env bash
# C++ test coverage for the standalone library (the doctest suite at tests/).
#
# Uses clang source-based instrumentation (-fprofile-instr-generate
# -fcoverage-mapping) plus llvm-cov / llvm-profdata. These ship with the clang
# toolchain, so there is nothing extra to install (no gcovr/lcov needed).
#
# Outputs:
#   * a per-file region/line/branch summary table on stdout
#   * an HTML report under build-cov/coverage-html/index.html
#
# Override the tool names if your clang is unsuffixed or a different version:
#   LLVM_COV=llvm-cov LLVM_PROFDATA=llvm-profdata dev/cpp-coverage.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LLVM_COV="${LLVM_COV:-llvm-cov-21}"
LLVM_PROFDATA="${LLVM_PROFDATA:-llvm-profdata-21}"

for tool in "$LLVM_COV" "$LLVM_PROFDATA"; do
  command -v "$tool" >/dev/null || {
    echo "error: '$tool' not found on PATH (set LLVM_COV / LLVM_PROFDATA)" >&2
    exit 127
  }
done

cmake --preset coverage
cmake --build --preset coverage --parallel

BIN="$(find build-cov -name misskappa_tests -type f -perm -u+x | head -1)"
[[ -n "$BIN" ]] || { echo "error: misskappa_tests not found under build-cov" >&2; exit 1; }

# Run the instrumented binary; one process -> one raw profile.
rm -f build-cov/cov.profraw build-cov/cov.profdata
LLVM_PROFILE_FILE="$ROOT/build-cov/cov.profraw" "$BIN"
"$LLVM_PROFDATA" merge -sparse build-cov/cov.profraw -o build-cov/cov.profdata

# Report only our own translation units: drop the test sources, the bundled
# doctest header, system headers, and Eigen.
IGNORE='(/tests/|third_party/|/usr/|[Ee]igen)'

echo
echo "================= C++ coverage (doctest suite) ================="
"$LLVM_COV" report "$BIN" \
  -instr-profile=build-cov/cov.profdata \
  -ignore-filename-regex="$IGNORE"

"$LLVM_COV" show "$BIN" \
  -instr-profile=build-cov/cov.profdata \
  -ignore-filename-regex="$IGNORE" \
  -format=html -output-dir=build-cov/coverage-html \
  -show-line-counts-or-regions >/dev/null

echo
echo "HTML report: build-cov/coverage-html/index.html"
