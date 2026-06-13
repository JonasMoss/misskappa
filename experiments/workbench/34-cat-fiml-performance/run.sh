#!/usr/bin/env bash
# Build and run the cat-FIML performance probe (three variants of
# src/estimate_fiml.cpp; see README.md). Eigen comes from the magmaan
# FetchContent checkout, matching build-opt's CMake cache.
set -euo pipefail
cd "$(dirname "$0")"

EIGEN="${EIGEN:-/home/jonas/Files/research/magmaan/build/fast/_deps/eigen3-src}"
FLAGS="-std=c++17 -O3 -g -DNDEBUG -DEIGEN_NO_EXCEPTIONS"
INC="-I../../../include -I../../../src -I. -I$EIGEN"

for v in prof fast sq; do
  g++ $FLAGS $INC bench.cpp "estimate_fiml_$v.cpp" -o "bench_$v"
done

REPS="${1:-2}"
for v in prof fast sq; do
  echo "=== bench_$v (reps=$REPS) ==="
  "./bench_$v" "$REPS"
done
