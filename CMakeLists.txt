cmake_minimum_required(VERSION 3.20)
project(poly_interpolation LANGUAGES CXX)

option(PI_ENABLE_WARNINGS "Enable warnings" ON)
option(PI_ENABLE_FLINT_BENCH "Build FLINT benchmark (requires libflint + libgmp)" ON)

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
  target_compile_options(interp_bench PRIVATE -O3 -DNDEBUG -march=native)
else()
  target_compile_options(interp_bench PRIVATE /O2)
  target_compile_definitions(interp_bench PRIVATE NDEBUG)
endif()

# ---- FLINT benchmark (external dependency) ----
if(PI_ENABLE_FLINT_BENCH)
  find_path(FLINT_INCLUDE_DIR flint/flint.h)
  find_library(FLINT_LIBRARY flint)
  find_library(GMP_LIBRARY gmp)
  if(NOT FLINT_INCLUDE_DIR OR NOT FLINT_LIBRARY OR NOT GMP_LIBRARY)
    message(WARNING "FLINT or GMP not found, skipping bench_flint target. Install libflint-dev/libgmp-dev or set PI_ENABLE_FLINT_BENCH=OFF to silence this.")
    set(PI_ENABLE_FLINT_BENCH OFF CACHE BOOL "" FORCE)
  endif()
endif()

if(PI_ENABLE_FLINT_BENCH)
  add_executable(bench_flint include/poly_interp/bench_flint.cpp)
  target_include_directories(bench_flint PRIVATE ${FLINT_INCLUDE_DIR})
  target_link_libraries(bench_flint PRIVATE poly_interp ${FLINT_LIBRARY} ${GMP_LIBRARY})

  if (NOT MSVC)
    target_compile_options(bench_flint PRIVATE -O3 -DNDEBUG -march=native)
  else()
    target_compile_options(bench_flint PRIVATE /O2)
    target_compile_definitions(bench_flint PRIVATE NDEBUG)
  endif()
endif()

# Convenience target: build + run bench to data/bench.csv
# Convenience targets: plot interp bench / plot flint bench (requires python + scripts/plot.py).
# Assumes corresponding CSV already exists.
add_custom_target(plot_interp
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/plots/interp
  COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot.py
          ${CMAKE_CURRENT_SOURCE_DIR}/data/bench.csv
          ${CMAKE_CURRENT_SOURCE_DIR}/plots/interp
          --format=interp
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)

add_custom_target(plot_flint
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/plots/flint
  COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot.py
          ${CMAKE_CURRENT_SOURCE_DIR}/data/bench_flint.csv
          ${CMAKE_CURRENT_SOURCE_DIR}/plots/flint
          --format=flint
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM
)
