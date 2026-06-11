# Tests are header-only doctest under -fno-exceptions.

add_library(doctest INTERFACE)
target_include_directories(doctest INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
target_compile_definitions(doctest INTERFACE
  DOCTEST_CONFIG_NO_EXCEPTIONS
  DOCTEST_CONFIG_NO_EXCEPTIONS_BUT_WITH_ALL_ASSERTS)

# Single test binary that links every unit/golden TU.
add_executable(misskappa_tests
  unit/tests_main.cpp
  unit/skeleton_test.cpp
  unit/loss_test.cpp
  unit/available_test.cpp
  unit/ipw_gwet_test.cpp
  unit/fiml_test.cpp
  unit/continuous_test.cpp
  unit/vector_test.cpp
  unit/normal_fiml_test.cpp
  unit/psd_inverse_test.cpp
  unit/counts_test.cpp
  unit/alpha_test.cpp
  unit/gwise_test.cpp
  unit/quadratic_test.cpp
  unit/fiml_counts_test.cpp
)
target_link_libraries(misskappa_tests PRIVATE misskappa doctest)
# psd_inverse_test.cpp white-box-tests a private detail_*.hpp header.
target_include_directories(misskappa_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME misskappa_tests COMMAND misskappa_tests)
