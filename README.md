# poly_interpolation

Lightweight C++17 header-only helpers for prime-field arithmetic and fast polynomial operations (NTT/CRT-backed multiplication, multipoint evaluation, and interpolation). Designed to handle large moduli (up to 64-bit primes such as Mersenne 2^61-1) and sizes up to 2^20 points.

## Requirements

- C++17 compiler with `unsigned __int128` support (GCC/Clang).
- CMake ≥ 3.20.
- (Optional, for FLINT baseline bench) `libflint-dev`, `libgmp-dev`.

## Building

```bash
cmake --preset release          # configure Ninja build in build/release
cmake --build --preset release --target interp_bench bench_flint
# or simply: make release
```

### Tests

```bash
cmake --build --preset release --target test
ctest --preset release
```

### Benchmarks

- Interp bench (project implementation): `make bench-interp` → writes `data/bench.csv`.
- Core ops bench (poly mul/div/mod, tree, batch inv, etc.): `make bench-ops` → writes `data/bench_ops.csv`.
- FLINT ops bench (nmod_poly versions of the same ops): `make bench-ops-flint` → writes `data/bench_ops_flint.csv` (requires FLINT/GMP).
- FLINT bench (requires FLINT/GMP): `make bench-flint` → writes `data/bench_flint.csv`.
- Plotting (requires python3 + pandas + matplotlib):
  - `make plot-interp` (reads `data/bench.csv`, outputs `plots/interp/`).
  - `make plot-ops` (reads `data/bench_ops.csv`, outputs `plots/ops/`, uses log-y).
  - `make plot-ops-flint` (reads `data/bench_ops_flint.csv`, outputs `plots/ops_flint/`, uses log-y).
  - `make plot-flint` (reads `data/bench_flint.csv`, outputs `plots/flint/`).
  - `python3 scripts/plot.py data/bench_ops_flint.csv plots/ops_flint --format=ops` (FLINT ops plots).

## Core APIs

All headers live under `include/poly_interp/`.

### `pf::FpCtx` and `pf::Fp` (`prime_field.hpp`)

- `FpCtx ctx(p)`: prime modulus context; provides `zero()`, `one()`, `from_uint()/from_int()`.
- Arithmetic: `add/sub/neg/mul/sqr/addmul/pow/inv/div`.
- Utilities: `batch_inv(std::vector<Fp>&)` for O(n) batch inversion.
- `Fp` is a thin wrapper holding `u64 v` (always canonicalized modulo `p`).

### `pf::FpPoly` (`fp_poly.hpp`)

Represents a univariate polynomial over a given `FpCtx`.

Key operations:

- Construction: `FpPoly(ctx)`, `FpPoly(ctx, std::vector<Fp>)`, `FpPoly(ctx, {u64...})`.
- Basic ops: `add/sub/mul` (NTT/CRT for large sizes), `scalar_mul`, `neg_poly`, `derivative`, `eval`.
- Division: `divrem`, `mod` with fast/slow paths.
- Multipoint evaluation: `multipoint_eval_naive(xs)`; `multipoint_eval(SubproductTree)`; convenience `multipoint_eval(xs)`.
- Interpolation: `interpolate_lagrange_naive`, `interpolate_subproduct_tree`.
- Helpers: `degree`, `is_zero`, `trim`, `coeff(i)`, `constant_term`, `leading_coeff`, `operator==/!=` and arithmetic operator overloads.

### `pf::FpPoly::SubproductTree`

- Built via `SubproductTree::build(FpCtx, xs)` to support fast multipoint evaluation and interpolation.
- Accessors: `root()`, `points`, `levels`, `n_points()`, `n_levels()`.
- Used by `FpPoly::multipoint_eval` and `FpPoly::interpolate_subproduct_tree`.

### NTT/CRT backend (internal)

- Five NTT primes (`kNTT`) with max_base ≥ 21 to cover 2^21-length convolutions; combined via specialized `CRT5Plan` (3+2 grouping).
- Convolution reuses buffers and thread-local scratch to minimize allocations.
  These are internal details; public APIs remain the simple `FpPoly::mul` interface.

## Usage Example

```cpp
#include <poly_interp/prime_field.hpp>
#include <poly_interp/fp_poly.hpp>
using namespace pf;

int main() {
    FpCtx F((1ULL<<61) - 1);           // M61
    FpPoly f(F, {1, 2, 3});            // 1 + 2x + 3x^2
    FpPoly g(F, {4, 5});               // 4 + 5x

    FpPoly h = f.mul(g);               // polynomial multiplication (NTT/CRT when large)
    Fp x = F.from_uint(7);
    Fp y = h.eval(x);                  // evaluate at x=7

    std::vector<Fp> xs = {F.from_uint(1), F.from_uint(2), F.from_uint(3)};
    auto tree = FpPoly::SubproductTree::build(F, xs);
    std::vector<Fp> ys = f.multipoint_eval(tree);         // f on xs
    FpPoly interp = FpPoly::interpolate_subproduct_tree(tree, ys); // recovers f
}
```

## Project Layout

- `include/poly_interp/prime_field.hpp` – prime-field context and element utilities.
- `include/poly_interp/fp_poly.hpp` – polynomial type and algorithms (NTT/CRT, multipoint eval, interpolation).
- `tests/` – unit tests for field operations, polynomials, and interpolation.
- `bench/` – benchmark driver; `data/` holds CSV output; `scripts/` contains plotting helper.
