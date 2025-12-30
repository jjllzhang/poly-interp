#include <poly_interp/prime_field.hpp>
#include <poly_interp/fp_poly.hpp>

#include <cassert>
#include <random>
#include <iostream>
#include <unordered_set>

using pf::Fp;
using pf::FpCtx;
using pf::FpPoly;
using pf::u64;

static pf::FpPoly mod_naive_ref(const pf::FpPoly& a, const pf::FpPoly& b) {
    const pf::FpCtx& F = a.ctx();
    // basic checks
    if (b.is_zero()) throw std::domain_error("mod_naive_ref: division by zero");
    if (a.is_zero()) return pf::FpPoly(F);

    const int degA = a.degree();
    const int degB = b.degree();
    if (degA < degB) return a;
    if (degB == 0) return pf::FpPoly(F); // mod constant -> 0

    // r = a (copy coeffs)
    std::vector<pf::Fp> r = a.coeffs();
    const pf::Fp lcB = b.leading_coeff();
    const bool monic = (lcB.v == 1);
    const pf::Fp inv_lcB = monic ? F.one() : F.inv(lcB);

    int degR = degA;

    auto get_b = [&](int i) -> pf::Fp {
        return b.coeff(static_cast<std::size_t>(i));
    };

    while (degR >= degB) {
        pf::Fp lead = r[static_cast<std::size_t>(degR)];
        if (lead.v != 0) {
            pf::Fp factor = monic ? lead : F.mul(lead, inv_lcB);
            int k = degR - degB;

            // r -= factor * x^k * b
            for (int i = 0; i < degB; ++i) {
                pf::Fp bi = get_b(i);
                if (bi.v == 0) continue;
                std::size_t idx = static_cast<std::size_t>(i + k);
                r[idx] = F.sub(r[idx], F.mul(factor, bi));
            }
            // cancel leading term
            r[static_cast<std::size_t>(degR)] = F.zero();
        } else {
            r[static_cast<std::size_t>(degR)] = F.zero();
        }

        --degR;
        while (degR >= 0 && r[static_cast<std::size_t>(degR)].v == 0) --degR;
    }

    if (degR < 0) return pf::FpPoly(F);
    r.resize(static_cast<std::size_t>(degR + 1));
    pf::FpPoly R(F, std::move(r));
    R.trim();
    return R;
}


static void test_trim_degree_basic() {
    FpCtx F(17);

    FpPoly a(F, {1, 2, 0, 0});
    assert(a.degree() == 1); // trim 后应该是 1

    FpPoly z(F);
    assert(z.degree() == -1);
    assert(z.is_zero());
}

static void test_add_sub_mul_derivative_eval() {
    FpCtx F(17);

    // f(x) = 3 + 5x + 0x^2 + 2x^3
    FpPoly f(F, {3, 5, 0, 2});
    assert(f.degree() == 3);

    // g(x) = 1 + x
    FpPoly g(F, {1, 1});

    // (1+x)^2 = 1 + 2x + x^2
    FpPoly gg = g.mul(g);
    assert(gg.degree() == 2);
    assert(gg.coeff(0).v == 1);
    assert(gg.coeff(1).v == 2);
    assert(gg.coeff(2).v == 1);

    // derivative: f'(x) = 5 + 6x^2 (mod 17)  因为 d(2x^3)=6x^2
    FpPoly df = f.derivative();
    assert(df.degree() == 2);
    assert(df.coeff(0).v == 5);
    assert(df.coeff(1).v == 0);
    assert(df.coeff(2).v == 6);

    // eval: f(2) = 3 + 5*2 + 2*8 = 3+10+16 = 29 = 12 mod 17
    Fp x = F.from_uint(2);
    Fp y = f.eval(x);
    assert(y.v == (29 % 17));
}

static void test_divrem_remainder_matches_eval() {
    FpCtx F(998244353ULL);

    // f(x) = x^3 + 2x + 1
    FpPoly f(F, {1, 2, 0, 1});

    // d(x) = x - 1
    FpPoly d(F);
    d.coeffs_mut().push_back(F.neg(F.from_uint(1)));
    d.coeffs_mut().push_back(F.one());
    d.trim();

    auto [q, r] = f.divrem(d);
    // r 应该是常数 f(1)=1+2+1=4
    assert(r.degree() <= 0);
    assert(r.constant_term().v == 4);
}

static void test_multipoint_eval_tree_vs_naive() {
    FpCtx F(998244353ULL);

    // random polynomial
    std::mt19937_64 rng(1234567);
    std::uniform_int_distribution<u64> dist(0, F.modulus() - 1);

    auto rand_fp = [&]() { return F.from_uint(dist(rng)); };

    for (int round = 0; round < 200; ++round) {
        // degree 0..50
        int deg = static_cast<int>(dist(rng) % 51);
        std::vector<Fp> coeffs(deg + 1);
        for (int i = 0; i <= deg; ++i) coeffs[i] = rand_fp();

        FpPoly f(F, coeffs);

        // points 1..32
        std::vector<Fp> xs;
        xs.reserve(32);
        for (int i = 0; i < 32; ++i) xs.push_back(rand_fp());

        auto ys_naive = f.multipoint_eval_naive(xs);
        auto ys_tree  = f.multipoint_eval(xs); // build tree internally

        assert(ys_naive.size() == ys_tree.size());
        for (std::size_t i = 0; i < ys_naive.size(); ++i) {
            assert(ys_naive[i] == ys_tree[i]);
        }
    }
}

static void test_interpolate_lagrange_naive() {
    FpCtx F(998244353ULL);
    std::mt19937_64 rng(20250101);
    std::uniform_int_distribution<u64> dist(0, F.modulus() - 1);

    auto rand_fp = [&]() { return F.from_uint(dist(rng)); };

    // n points
    const std::size_t n = 16;

    // generate distinct x_i (非常重要：否则 denom=0)
    std::vector<Fp> xs;
    xs.reserve(n);
    std::unordered_set<u64> used;
    while (xs.size() < n) {
        Fp x = rand_fp();
        if (used.insert(x.v).second) xs.push_back(x);
    }

    // pick a random polynomial deg < n
    const int deg = 7;
    std::vector<Fp> coeffs(deg + 1);
    for (int i = 0; i <= deg; ++i) coeffs[i] = rand_fp();
    FpPoly f(F, coeffs);

    // ys = f(xs)
    std::vector<Fp> ys = f.multipoint_eval_naive(xs);

    // interpolate
    FpPoly g = FpPoly::interpolate_lagrange_naive(F, xs, ys);

    // verify on the interpolation points
    auto ys2 = g.multipoint_eval_naive(xs);
    for (std::size_t i = 0; i < n; ++i) {
        assert(ys2[i] == ys[i]);
    }

    // verify on some extra random points
    for (int t = 0; t < 50; ++t) {
        Fp x = rand_fp();
        assert(f.eval(x) == g.eval(x));
    }
}

static void test_interpolate_subproduct_tree_vs_naive() {
    FpCtx F(998244353ULL);
    std::mt19937_64 rng(20251229);
    std::uniform_int_distribution<u64> dist(0, F.modulus() - 1);

    auto rand_fp = [&]() { return F.from_uint(dist(rng)); };

    for (int round = 0; round < 200; ++round) {
        const std::size_t n = 32;

        // 生成互异点 xs
        std::vector<Fp> xs;
        xs.reserve(n);
        std::unordered_set<u64> used;
        while (xs.size() < n) {
            Fp x = rand_fp();
            if (used.insert(x.v).second) xs.push_back(x);
        }

        // 随机多项式 deg < n
        int deg = 10;
        std::vector<Fp> coeffs(deg + 1);
        for (int i = 0; i <= deg; ++i) coeffs[i] = rand_fp();
        FpPoly f(F, coeffs);

        // 计算 ys
        std::vector<Fp> ys = f.multipoint_eval_naive(xs);

        // 构建树
        auto tree = FpPoly::SubproductTree::build(F, xs);

        // 子乘积树插值
        FpPoly g = FpPoly::interpolate_subproduct_tree(tree, ys);

        // 对拍：在所有插值点上应该一致
        auto ys2 = g.multipoint_eval(tree);
        for (std::size_t i = 0; i < n; ++i) {
            assert(ys2[i] == ys[i]);
        }

        // 再随机测几个点
        for (int t = 0; t < 20; ++t) {
            Fp x = rand_fp();
            assert(f.eval(x) == g.eval(x));
        }
    }
}

static void test_interpolate_roundtrip(u64 p) {
    pf::FpCtx F(p);
    std::mt19937_64 rng(777);

    const std::size_t n = 4096;
    std::vector<pf::Fp> xs(n), ys(n);

    for (std::size_t i = 0; i < n; ++i) xs[i] = F.from_uint((pf::u64)(i + 1));
    for (std::size_t i = 0; i < n; ++i) ys[i] = F.from_uint((pf::u64)rng());

    auto tree = pf::FpPoly::SubproductTree::build(F, xs);
    pf::FpPoly poly = pf::FpPoly::interpolate_subproduct_tree(tree, ys);

    std::vector<pf::Fp> got = poly.multipoint_eval(tree);
    assert(got.size() == ys.size());
    for (std::size_t i = 0; i < n; ++i) assert(got[i] == ys[i]);
}

static void test_poly_operator_overloads() {
    FpCtx F(17);
    FpPoly f(F, {3, 5, 0, 2});
    FpPoly g(F, {1, 1});

    assert((f + g) == f.add(g));
    assert((f - g) == f.sub(g));
    assert((f * g) == f.mul(g));
    assert((-f) == f.neg_poly());

    FpPoly h = f;
    h += g;
    assert(h == f.add(g));

    h = f;
    h -= g;
    assert(h == f.sub(g));

    h = f;
    h *= g;
    assert(h == f.mul(g));

    // scalar
    pf::Fp k = F.from_uint(7);
    assert((f * k) == f.scalar_mul(k));
    assert((k * f) == f.scalar_mul(k));

    // operator()
    pf::Fp x = F.from_uint(2);
    assert(f(x) == f.eval(x));
}

static void test_remainder_tree_leaf_values() {
    FpCtx F(998244353ULL);

    // f(x)=1+2x+3x^2
    FpPoly f(F, {1, 2, 3});

    std::vector<pf::Fp> xs;
    for (u64 i = 1; i <= 32; ++i) xs.push_back(F.from_uint(i));

    auto tree = FpPoly::SubproductTree::build(F, xs);
    auto rem  = f.remainder_tree(tree);

    assert(!rem.empty());
    assert(rem[0].size() == xs.size());

    for (std::size_t i = 0; i < xs.size(); ++i) {
        pf::Fp yi = f.eval(xs[i]);
        pf::Fp ri = rem[0][i].constant_term();
        assert(yi == ri);
    }
}

static void test_batch_inv() {
    FpCtx F(998244353ULL);
    std::vector<Fp> v;
    for (u64 i = 1; i <= 1000; ++i) v.push_back(F.from_uint(i));

    std::vector<Fp> inv1 = v;
    F.batch_inv(inv1);

    for (std::size_t i = 0; i < v.size(); ++i) {
        assert(F.mul(v[i], inv1[i]) == F.one());
    }
}

static pf::FpPoly mul_naive_ref(const pf::FpPoly& a, const pf::FpPoly& b) {
    const pf::FpCtx& F = a.ctx();
    if (a.is_zero() || b.is_zero()) return pf::FpPoly(F);
    const std::size_t n = a.coeffs().size();
    const std::size_t m = b.coeffs().size();

    std::vector<pf::Fp> out(n + m - 1, F.zero());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            out[i + j] = F.add(out[i + j], F.mul(a.coeffs()[i], b.coeffs()[j]));
        }
    }

    pf::FpPoly r(F, std::move(out));
    r.trim();
    return r;
}


static void test_mul_ntt_crt_correctness() {
    pf::FpCtx F((1ULL<<61) - 1); // M61
    std::mt19937_64 rng(123);

    for (int round = 0; round < 50; ++round) {
        int n = 800 + (int)(rng() % 400);  // 让它超过阈值，触发 NTT/CRT
        int m = 700 + (int)(rng() % 400);

        std::vector<pf::Fp> A(n), B(m);
        for (int i = 0; i < n; ++i) A[i] = F.from_uint((pf::u64)rng());
        for (int j = 0; j < m; ++j) B[j] = F.from_uint((pf::u64)rng());

        pf::FpPoly a(F, A), b(F, B);
        auto c1 = a.mul(b);            // 你的实现（可能走 NTT/CRT）
        auto c2 = mul_naive_ref(a, b); // 参考朴素

        assert(c1 == c2);
    }
}

static void test_mul_eval_consistency(u64 p) {
    pf::FpCtx F(p);
    std::mt19937_64 rng(123456);

    auto rndFp = [&]() -> pf::Fp { return F.from_uint((pf::u64)rng()); };

    for (int tc = 0; tc < 80; ++tc) {
        // 让规模跨过阈值，逼迫走 NTT/CRT 路径
        std::size_t n = 400 + (rng() % 800);
        std::size_t m = 400 + (rng() % 800);

        std::vector<pf::Fp> A(n), B(m);
        for (std::size_t i = 0; i < n; ++i) A[i] = rndFp();
        for (std::size_t j = 0; j < m; ++j) B[j] = rndFp();

        // 偶尔用极端系数 (p-1) 压一压范围
        if (tc % 10 == 0) {
            for (auto& x : A) x = F.sub(F.zero(), F.one());
            for (auto& x : B) x = F.sub(F.zero(), F.one());
        }

        pf::FpPoly a(F, A), b(F, B);
        pf::FpPoly c = a.mul(b);

        // 多取几个随机点
        for (int rep = 0; rep < 16; ++rep) {
            pf::Fp x = F.from_uint((pf::u64)(rng() % 100000 + 1));
            pf::Fp lhs = c.eval(x);
            pf::Fp rhs = F.mul(a.eval(x), b.eval(x));
            assert(lhs == rhs);
        }
    }
}

static void test_fast_mod_matches_slow() {
    pf::FpCtx F((1ULL<<61)-1);
    std::mt19937_64 rng(12345);
    auto rnd = [&]() { return F.from_uint((pf::u64)rng()); };

    for (int round = 0; round < 20; ++round) {
        int n = 2000 + (int)(rng() % 500);  // 让它触发 fast
        int m = 800  + (int)(rng() % 200);

        std::vector<pf::Fp> A(n), B(m);
        for (int i = 0; i < n; ++i) A[i] = rnd();
        for (int j = 0; j < m; ++j) B[j] = rnd();

        // 保证 divisor 首项非零（trim 后的 leading coeff）
        B.back() = F.one();

        pf::FpPoly a(F, A), b(F, B);

        auto r_fast = a.mod(b);        // dispatch -> fast
        auto r_slow = mod_naive_ref(a, b);

        assert(r_fast == r_slow);
    }
}




int main() {
    test_trim_degree_basic();
    test_add_sub_mul_derivative_eval();
    test_divrem_remainder_matches_eval();
    test_multipoint_eval_tree_vs_naive();
    test_interpolate_lagrange_naive();
    test_interpolate_subproduct_tree_vs_naive();
    test_interpolate_roundtrip((1ULL<<31) - 1);
    test_interpolate_roundtrip((1ULL<<61) - 1);
    test_poly_operator_overloads();
    test_remainder_tree_leaf_values();
    test_batch_inv();
    test_mul_eval_consistency((1ULL<<31) - 1);
    test_mul_eval_consistency((1ULL<<61) - 1);
    test_mul_ntt_crt_correctness();
    test_fast_mod_matches_slow();

    std::cout << "All fp_poly tests passed.\n";
    return 0;
}
