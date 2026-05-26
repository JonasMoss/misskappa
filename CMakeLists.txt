cmake_minimum_required(VERSION 3.20)
project(misskappa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Non-negotiables: no exceptions, no RTTI. Eigen runs under EIGEN_NO_EXCEPTIONS.
add_compile_options(-fno-exceptions -fno-rtti -Wall -Wextra -Wpedantic)
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
  src/estimate_raw.cpp
)
# Allow private detail_*.hpp to be included with relative paths from src/.
target_include_directories(misskappa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(misskappa PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(misskappa PUBLIC Eigen3::Eigen)

# Tests
enable_testing()
add_subdirectory(tests)
