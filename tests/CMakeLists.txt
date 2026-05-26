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
)
target_link_libraries(misskappa_tests PRIVATE misskappa doctest)
add_test(NAME misskappa_tests COMMAND misskappa_tests)
