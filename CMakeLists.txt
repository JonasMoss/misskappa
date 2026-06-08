cmake_minimum_required(VERSION 3.20)
project(misskappa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -fPIC so libmisskappa.a can link into the R package shared object.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Non-negotiables: no exceptions, no RTTI. Eigen runs under EIGEN_NO_EXCEPTIONS.
# (-fno-exceptions / -fno-rtti are applied only to the library target below, so
# the R package's Rcpp glue can compile with the defaults R expects.)
add_compile_options(-Wall -Wextra -Wpedantic)
add_compile_definitions(EIGEN_NO_EXCEPTIONS)

# Eigen: prefer system / submodule. Fail fast if not found rather than
# silently fetching, so the toolchain stays explicit.
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
    HINTS /usr/include/eigen3 /usr/local/include/eigen3
    REQUIRED)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
endif()

add_library(misskappa STATIC
  src/loss.cpp
  src/detail_inverse_weights.cpp
  src/estimate_raw.cpp
  src/estimate_continuous.cpp
  src/estimate_vector.cpp
  src/estimate_counts.cpp
  src/estimate_alpha.cpp
  src/estimate_gwise.cpp
  src/estimate_quadratic.cpp
  src/estimate_fiml.cpp
  src/estimate_fiml_counts.cpp
)
# Allow private detail_*.hpp to be included with relative paths from src/.
target_include_directories(misskappa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(misskappa PRIVATE -fno-exceptions -fno-rtti)
target_include_directories(misskappa PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(misskappa PUBLIC Eigen3::Eigen)

# Tests
enable_testing()
add_subdirectory(tests)
