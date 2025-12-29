#include <poly_interp/prime_field.hpp>
#include <poly_interp/fp_poly.hpp>

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using pf::Fp;
using pf::FpCtx;
using pf::FpPoly;
using pf::u64;

static u64 mersenne31() { return (1ULL << 31) - 1; }
static u64 mersenne61() { return (1ULL << 61) - 1; }

static std::vector<Fp> make_distinct_xs_pow2(const FpCtx& F, int k) {
    // x_i = i, i=0..n-1  (保证互异)
    const std::size_t n = (std::size_t)1ULL << (unsigned)k;
    std::vector<Fp> xs(n);
    for (std::size_t i = 0; i < n; ++i) xs[i] = F.from_uint((u64)i);
    return xs;
}

static std::vector<Fp> random_poly_coeffs(const FpCtx& F, std::size_t deg, std::mt19937_64& rng) {
    std::uniform_int_distribution<u64> dist(0, F.modulus() - 1);
    std::vector<Fp> c(deg + 1);
    for (std::size_t i = 0; i <= deg; ++i) c[i] = F.from_uint(dist(rng));
    // 保证最高项非 0（让 degree 固定）
    if (deg > 0 && c[deg].v == 0) c[deg] = F.one();
    return c;
}

static void verify_on_points(const FpPoly& f,
                             const std::vector<Fp>& xs,
                             const std::vector<Fp>& ys,
                             std::size_t samples,
                             std::mt19937_64& rng) {
    const std::size_t n = xs.size();
    assert(ys.size() == n);

    if (samples >= n) {
        // 全量验证
        for (std::size_t i = 0; i < n; ++i) {
            Fp got = f.eval(xs[i]);
            if (got != ys[i]) {
                std::cerr << "Mismatch at i=" << i << " got=" << got.v << " expected=" << ys[i].v << "\n";
                std::abort();
            }
        }
        return;
    }

    // 抽样验证
    std::uniform_int_distribution<std::size_t> pick(0, n - 1);
    for (std::size_t t = 0; t < samples; ++t) {
        std::size_t i = pick(rng);
        Fp got = f.eval(xs[i]);
        if (got != ys[i]) {
            std::cerr << "Mismatch at sampled i=" << i << " got=" << got.v << " expected=" << ys[i].v << "\n";
            std::abort();
        }
    }
}

static void test_interp_roundtrip_poly(const FpCtx& F,
                                       int k,
                                       std::size_t deg,
                                       std::size_t verify_samples,
                                       u64 seed) {
    std::mt19937_64 rng(seed ^ (F.modulus() + 0x9e3779b97f4a7c15ULL) ^ (u64)k);

    const auto xs = make_distinct_xs_pow2(F, k);
    const std::size_t n = xs.size();
    // 生成一个随机多项式 g（degree <= deg）
    auto coeffs = random_poly_coeffs(F, deg, rng);
    FpPoly g(F, coeffs);

    // ys = g(xs)  (用 multipoint_eval/tree 或者逐点 eval 都可以；逐点更“独立”)
    std::vector<Fp> ys(n, F.zero());
    for (std::size_t i = 0; i < n; ++i) ys[i] = g.eval(xs[i]);

    // 插值：f 应该在这些点上等于 g
    auto tree = FpPoly::SubproductTree::build(F, xs);
    FpPoly f = FpPoly::interpolate_subproduct_tree(tree, ys);

    // 验证：点上相等（必要且充分）
    verify_on_points(f, xs, ys, verify_samples, rng);

    // 额外：如果 deg < n，理论上 f == g (作为多项式) 也应成立
    // 但考虑你的内部 trim / 表示差异，这里用 equals 即可
    // 若你希望更严，可以抽样比较更多随机点：
    if (deg + 1 <= n) {
        // 随机再测一些点（不在 xs 集合里）
        std::uniform_int_distribution<u64> dist(0, F.modulus() - 1);
        for (int t = 0; t < 32; ++t) {
            Fp x = F.from_uint(dist(rng));
            if (f.eval(x) != g.eval(x)) {
                std::cerr << "Mismatch at extra random x, got=" << f.eval(x).v
                          << " expected=" << g.eval(x).v << "\n";
                std::abort();
            }
        }
    }
}

static void test_interp_random_ys(const FpCtx& F,
                                 int k,
                                 std::size_t verify_samples,
                                 u64 seed) {
    std::mt19937_64 rng(seed ^ (F.modulus() * 1315423911ULL) ^ (u64)k);
    const auto xs = make_distinct_xs_pow2(F, k);
    const std::size_t n = xs.size();

    // 随机 ys（更贴近你的 bench）
    std::uniform_int_distribution<u64> dist(0, F.modulus() - 1);
    std::vector<Fp> ys(n);
    for (std::size_t i = 0; i < n; ++i) ys[i] = F.from_uint(dist(rng));

    auto tree = FpPoly::SubproductTree::build(F, xs);
    FpPoly f = FpPoly::interpolate_subproduct_tree(tree, ys);

    verify_on_points(f, xs, ys, verify_samples, rng);
}

int main() {
    // 覆盖 M31/M61
    const u64 seeds[] = {1, 2, 3, 12345, 7777777};

    for (u64 p : {mersenne31(), mersenne61()}) {
        FpCtx F(p);
        std::cerr << "[interp_correctness] prime=" << p << "\n";

        // 小规模：全量验证
        for (int k : {8, 9, 10, 11}) { // n=256..2048
            for (u64 sd : seeds) {
                // roundtrip：随机多项式 -> 点值 -> 插值 -> 验证
                test_interp_roundtrip_poly(F, k, /*deg=*/std::min<std::size_t>((1u<<k) - 1, 256), /*verify_samples=*/(1u<<k), sd);
                // 随机 ys 也测一下
                test_interp_random_ys(F, k, /*verify_samples=*/(1u<<k), sd ^ 0xabcdef);
            }
        }

        // 中/大规模：抽样验证（避免测试跑太久）
        for (int k : {12, 13, 14, 15, 16}) { // n=4096..65536
            for (u64 sd : seeds) {
                test_interp_random_ys(F, k, /*verify_samples=*/128, sd);
                test_interp_roundtrip_poly(F, k, /*deg=*/256, /*verify_samples=*/128, sd ^ 0x12345678);
            }
        }
    }

    std::cerr << "All interp_correctness tests passed.\n";
    return 0;
}
