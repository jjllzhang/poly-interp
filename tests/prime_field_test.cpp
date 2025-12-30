#include <poly_interp/prime_field.hpp>

#include <cassert>
#include <random>
#include <iostream>
#include <limits>

using pf::Fp;
using pf::FpCtx;
using pf::u64;
using pf::u128;

static void test_small_handcrafted() {
    FpCtx F(17);

    Fp a = F.from_uint(5);
    Fp b = F.from_uint(13);

    // 5 + 13 = 18 = 1 (mod 17)
    assert(F.add(a, b).v == 1);

    // 5 - 13 = -8 = 9 (mod 17)
    assert(F.sub(a, b).v == 9);

    // 5 * 13 = 65 = 14 (mod 17)
    assert(F.mul(a, b).v == 14);

    // inv(5) = 7 because 5*7=35=1 mod 17
    assert(F.mul(a, F.inv(a)).v == 1);
    assert(F.inv(a).v == 7);

    // pow
    assert(F.pow(a, 0).v == 1);
    assert(F.pow(a, 1).v == a.v);
}

static void test_mod2() {
    FpCtx F(2);
    Fp zero = F.zero();
    Fp one  = F.one();

    assert(one.v == 1);
    assert(F.add(one, one).v == 0);
    assert(F.mul(one, one).v == 1);
    assert(F.inv(one).v == 1);

    // inv(0) should throw
    bool threw = false;
    try { (void)F.inv(zero); }
    catch (const std::domain_error&) { threw = true; }
    assert(threw);
}

static void test_random_identities(u64 p, int rounds, std::uint64_t seed) {
    FpCtx F(p);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<u64> dist(0, p - 1);

    auto rand_elem = [&]() -> Fp { return F.from_uint(dist(rng)); };
    auto rand_nonzero = [&]() -> Fp {
        while (true) {
            Fp x = rand_elem();
            if (!F.is_zero(x)) return x;
        }
    };

    for (int i = 0; i < rounds; ++i) {
        Fp a = rand_elem();
        Fp b = rand_elem();
        Fp c = rand_elem();

        // 结果必须在 [0, p)
        auto chk = [&](Fp x) { assert(x.v < p); };

        chk(a); chk(b); chk(c);

        // 加法交换律/结合律
        chk(F.add(a, b));
        assert(F.add(a, b) == F.add(b, a));
        assert(F.add(F.add(a, b), c) == F.add(a, F.add(b, c)));

        // 加法单位元/逆元
        assert(F.add(a, F.zero()) == a);
        assert(F.add(a, F.neg(a)) == F.zero());

        // 乘法交换律/结合律
        assert(F.mul(a, b) == F.mul(b, a));
        assert(F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c)));

        // 乘法单位元
        assert(F.mul(a, F.one()) == a);

        // 分配律：a*(b+c)=a*b+a*c
        assert(F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c)));

        // 减法一致性：(a-b)+b=a
        assert(F.add(F.sub(a, b), b) == a);

        // 非零元素逆元：a*inv(a)=1
        Fp x = rand_nonzero();
        assert(F.mul(x, F.inv(x)) == F.one());

        // 除法一致性：(a/b)*b = a (b!=0)
        Fp y = rand_nonzero();
        assert(F.mul(F.div(a, y), y) == a);

        // Fermat: x^(p-1)=1 for x!=0 (可选但很好用)
        // 注意 p=2 时 p-1=1 仍成立
        assert(F.pow(x, p - 1) == F.one());
    }
}

static void test_mersenne_mul(u64 p, int rounds, std::uint64_t seed) {
    FpCtx F(p);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<u64> dist(0, p - 1);

    auto mul_ref = [&](u64 a, u64 b) -> u64 {
        return static_cast<u64>(static_cast<u128>(a) * static_cast<u128>(b) % static_cast<u128>(p));
    };

    for (int i = 0; i < rounds; ++i) {
        u64 a = dist(rng);
        u64 b = dist(rng);
        u64 ref = mul_ref(a, b);
        Fp got = F.mul(F.from_uint(a), F.from_uint(b));
        if (got.v != ref) {
            std::cerr << "Mersenne mul mismatch: p=" << p
                      << ", a=" << a << ", b=" << b
                      << ", got=" << got.v << ", expected=" << ref << "\n";
            std::abort();
        }
    }

    // from_uint reduction for values >= p
    for (int i = 0; i < 1000; ++i) {
        u64 raw = dist(rng) + p; // in [p, 2p)
        u64 ref = raw % p;
        if (F.from_uint(raw).v != ref) {
            std::cerr << "from_uint mismatch for Mersenne p=" << p
                      << ", raw=" << raw << ", expected=" << ref << "\n";
            std::abort();
        }
    }
}

static void test_mersenne_primes() {
    const u64 m31 = (1ULL << 31) - 1;
    const u64 m61 = (1ULL << 61) - 1;

    test_mersenne_mul(m31, 20000, 0x12345678ULL);
    test_mersenne_mul(m61, 20000, 0x9abcdef0ULL);
}

int main() {
    test_small_handcrafted();
    test_mod2();

    // 常见素数：998244353（NTT 常用）
    test_random_identities(998244353ULL, 20000, 1234567);

    // 也可以测一个较大的素数（2^61-1 是著名梅森素数）
    test_random_identities(2305843009213693951ULL, 20000, 7654321);

    // 梅森素数特化路径
    test_mersenne_primes();

    std::cout << "All prime field tests passed.\n";
    return 0;
}
