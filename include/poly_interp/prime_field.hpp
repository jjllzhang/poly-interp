#ifndef PRIME_FIELD_HPP
#define PRIME_FIELD_HPP

#include <cstdint>
#include <stdexcept>
#include <ostream>
#include <vector>
#include <cstddef>


namespace pf {

// 说明：这里用 unsigned __int128 来保证 64-bit 模数下的正确性。
// 如果你用 MSVC，需要改成 boost::multiprecision 或自写 128 位乘法/取模。
#if !defined(__SIZEOF_INT128__)
#  error "prime_field.hpp requires compiler support for unsigned __int128 (GCC/Clang)."
#endif

using u64  = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

// 素域元素：仅存值 v，约定始终保持 0 <= v < p（规范表示）。
struct Fp {
    u64 v = 0;
};

// 仅做调试输出（可选）。
inline std::ostream& operator<<(std::ostream& os, const Fp& x) {
    return os << x.v;
}

// 素域上下文：只保存模数 p（假设 p 是素数）。
struct FpCtx {
    u64 p;
    bool is_mersenne_ = false;
    unsigned mersenne_k_ = 0; // p = 2^k - 1
    u64 mersenne_mask_ = 0;

    // 不做素性测试；只做最基本合法性检查。
    explicit FpCtx(u64 prime_mod) : p(prime_mod) {
        if (p < 2) {
            throw std::invalid_argument("FpCtx: modulus p must be >= 2");
        }

        // Detect p = 2^k - 1 (Mersenne); used for fast reduction paths.
        is_mersenne_ = ((p & (p + 1)) == 0);
        if (is_mersenne_) {
            mersenne_mask_ = p; // == (1ULL << k) - 1
            mersenne_k_ = 64U - static_cast<unsigned>(__builtin_clzll(p));
            if (mersenne_k_ >= 64U) { // avoid undefined shift when p == 2^64-1
                is_mersenne_ = false;
                mersenne_k_ = 0;
                mersenne_mask_ = 0;
            }
        }
    }

    u64 modulus() const noexcept { return p; }

    // 常量
    Fp zero() const noexcept { return Fp{0}; }
    Fp one()  const noexcept { return Fp{1 % p}; }

    // 将整数规约到 [0, p)
    Fp from_uint(u64 x) const noexcept {
        if (is_mersenne_) {
            return Fp{ reduce_mersenne(static_cast<u128>(x)) };
        }
        return Fp{ x % p };
    }

    // 将有符号整数规约到 [0, p)
    Fp from_int(std::int64_t x) const noexcept {
        if (is_mersenne_) {
            i128 xi = static_cast<i128>(x);
            bool neg = xi < 0;
            u128 mag = neg ? static_cast<u128>(-xi) : static_cast<u128>(xi);
            u64 r = reduce_mersenne(mag);
            if (neg && r != 0) r = p - r;
            return Fp{r};
        }

        i128 r = static_cast<i128>(x) % static_cast<i128>(p);
        if (r < 0) r += static_cast<i128>(p);
        return Fp{ static_cast<u64>(r) };
    }

    // 基本谓词
    bool is_zero(Fp a) const noexcept { return a.v == 0; }
    bool eq(Fp a, Fp b) const noexcept { return a.v == b.v; }

    // Optimized modular addition for canonical operands (0 <= a,b < p).
    // Avoids division/mod and is overflow-safe for any uint64 p.
    Fp add(Fp a, Fp b) const noexcept {
        // r = a + b mod p
        // If a >= p-b then a+b-p else a+b
        const u64 threshold = p - b.v; // in [1..p]
        const u64 r = (a.v >= threshold) ? (a.v - threshold) : (a.v + b.v);
        return Fp{r};
    }

    // Optimized modular subtraction for canonical operands (0 <= a,b < p).
    Fp sub(Fp a, Fp b) const noexcept {
        if (a.v >= b.v) return Fp{a.v - b.v};
        // Here a < b, so a + (p-b) < p, no overflow.
        return Fp{ a.v + (p - b.v) };
    }


    // -a (mod p)
    Fp neg(Fp a) const noexcept {
        if (a.v == 0) return a;
        return Fp{ p - a.v };
    }

    // a * b (mod p)
    Fp mul(Fp a, Fp b) const noexcept {
        u128 t = static_cast<u128>(a.v) * static_cast<u128>(b.v);
        u64 r = is_mersenne_ ? reduce_mersenne(t) : static_cast<u64>(t % p);
        return Fp{r};
    }

    // a^2 (mod p)
    Fp sqr(Fp a) const noexcept {
        return mul(a, a);
    }

    // r + a*b (mod p) ——多项式/向量运算里很常用（先给个正确版本）
    Fp addmul(Fp r, Fp a, Fp b) const noexcept {
        return add(r, mul(a, b));
    }

    // 快速幂：a^e (mod p)
    Fp pow(Fp a, u64 e) const noexcept {
        Fp base = a;
        Fp res  = one();
        u64 exp = e;

        while (exp > 0) {
            if (exp & 1ULL) res = mul(res, base);
            exp >>= 1ULL;
            if (exp) base = sqr(base);
        }
        return res;
    }

    // 逆元：a^{-1} (mod p)
    // 由于 p 是素数，使用费马小定理：a^(p-2) mod p
    Fp inv(Fp a) const {
        if (a.v == 0) {
            throw std::domain_error("FpCtx::inv: inverse of zero");
        }
        // p=2 时，p-2=0，inv(1)=1 也成立
        return pow(a, p - 2);
    }

    // 除法：a / b = a * inv(b)
    Fp div(Fp a, Fp b) const {
        return mul(a, inv(b));
    }

    // In-place batch inversion:
    // input:  vec[i] in F_p, all must be non-zero
    // output: vec[i] = inv(vec[i])
    //
    // Complexity: O(n) mul + 1 inv
    void batch_inv(std::vector<Fp>& vec) const {
        const std::size_t n = vec.size();
        if (n == 0) return;

        // Check non-zero and canonicalize (defensive)
        for (auto& x : vec) {
            x.v %= p;
            if (x.v == 0) {
                throw std::domain_error("FpCtx::batch_inv: zero element in batch");
            }
        }

        // prefix products: pref[i] = vec[0]*...*vec[i]
        std::vector<Fp> pref(n);
        pref[0] = vec[0];
        for (std::size_t i = 1; i < n; ++i) {
            pref[i] = mul(pref[i - 1], vec[i]);
        }

        // inv_total = 1 / (vec[0]*...*vec[n-1])
        Fp inv_total = inv(pref[n - 1]);

        // Backward pass:
        // vec[i] = inv_total * pref[i-1]
        // inv_total *= old vec[i]
        for (std::size_t i = n; i-- > 0;) {
            Fp before = (i == 0) ? one() : pref[i - 1];
            Fp old = vec[i];

            vec[i] = mul(inv_total, before);
            inv_total = mul(inv_total, old);
        }
    }

private:
    // Fast reduction for Mersenne prime p = 2^k - 1.
    inline u64 reduce_mersenne(u128 x) const noexcept {
        // Fold: (low k bits) + (rest), repeat once; final result < 2^k + 1 <= 2p
        const u128 mask = static_cast<u128>(mersenne_mask_);
        const unsigned k = mersenne_k_;

        u128 t = (x & mask) + (x >> k);
        u64 r = static_cast<u64>((t & mask) + (t >> k));
        if (r >= mersenne_mask_) r -= mersenne_mask_;
        return r;
    }
};

// 方便用的比较运算（不依赖 ctx；仅比较存储值）
inline bool operator==(Fp a, Fp b) noexcept { return a.v == b.v; }
inline bool operator!=(Fp a, Fp b) noexcept { return a.v != b.v; }

} // namespace pf

#endif // PRIME_FIELD_HPP
