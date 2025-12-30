cmake_minimum_required(VERSION 3.20)
project(poly_interpolation LANGUAGES CXX)

option(PI_ENABLE_WARNINGS "Enable warnings" ON)

# ---- Header-only library target ----
add_library(poly_interp INTERFACE)
target_include_directories(poly_interp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(poly_interp INTERFACE cxx_std_17)


if(PI_ENABLE_WARNINGS)
  if (MSVC)
    target_compile_options(poly_interp INTERFACE /W4)
  else()
    target_compile_options(poly_interp INTERFACE -Wall -Wextra -Wpedantic -Wno-pedantic)
  endif()
endif()


# ---- Tests ----
include(CTest)

add_executable(prime_field_test tests/prime_field_test.cpp)
target_link_libraries(prime_field_test PRIVATE poly_interp)
add_test(NAME prime_field_test COMMAND prime_field_test)

add_executable(fp_poly_test tests/fp_poly_test.cpp)
target_link_libraries(fp_poly_test PRIVATE poly_interp)
add_test(NAME fp_poly_test COMMAND fp_poly_test)

add_executable(interp_correctness_test tests/interp_correctness_test.cpp)
target_link_libraries(interp_correctness_test PRIVATE poly_interp)
add_test(NAME interp_correctness_test COMMAND interp_correctness_test)

if (NOT MSVC)
  target_compile_options(prime_field_test PRIVATE -Wno-unused-but-set-variable -Wno-unused-parameter)
  target_compile_options(fp_poly_test PRIVATE -Wno-unused-but-set-variable -Wno-unused-parameter)
  target_compile_options(interp_correctness_test PRIVATE -Wno-unused-but-set-variable -Wno-unused-parameter)
endif()


# ---- Benchmark ----
add_executable(interp_bench bench/interp_bench.cpp)
target_link_libraries(interp_bench PRIVATE poly_interp)

# Bench: prefer Release-like flags
if (NOT MSVC)
  target_compile_options(interp_bench PRIVATE -O3 -DNDEBUG)
else()
  target_compile_options(interp_bench PRIVATE /O2)
  target_compile_definitions(interp_bench PRIVATE NDEBUG)
endif()

# Convenience target: build + run bench to data/bench.csv
add_custom_target(run_bench
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/data
  COMMAND interp_bench --prime=all --min_pow=10 --max_pow=20 --repeats=3 --csv=1
          > ${CMAKE_CURRENT_SOURCE_DIR}/data/bench.csv
  DEPENDS interp_bench
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
)

# Convenience target: plot (requires python + scripts/plot.py)
add_custom_target(plot
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/plots
  COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot.py
          ${CMAKE_CURRENT_SOURCE_DIR}/data/bench.csv
          ${CMAKE_CURRENT_SOURCE_DIR}/plots
  DEPENDS run_bench
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
