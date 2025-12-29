#ifndef FP_POLY_HPP
#define FP_POLY_HPP

#include <poly_interp/prime_field.hpp>


#include <vector>
#include <initializer_list>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <utility>
#include <ostream>
#include <array>
#include <cassert>
#include <cstdint>


namespace pf {

namespace detail {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct NTTPrime {
    u32 mod;
    u32 primitive_root;
    int max_base; // supports length up to 2^max_base
};

// 6 primes => product ~ 2^178, 足够覆盖 n<=2^20, p<2^64 的系数上界（保守）
static constexpr std::array<NTTPrime, 6> kNTT = {{
    {167772161u, 3u, 25},     // 5*2^25+1
    {469762049u, 3u, 26},     // 7*2^26+1
    {754974721u, 11u, 24},    // 45*2^24+1
    {998244353u, 3u, 23},     // 119*2^23+1
    {1224736769u, 3u, 24},    // 73*2^24+1
    {1004535809u, 3u, 21},    // 479*2^21+1
}};

inline u32 add_mod(u32 a, u32 b, u32 mod) {
    u32 s = a + b;
    if (s >= mod) s -= mod;
    return s;
}
inline u32 sub_mod(u32 a, u32 b, u32 mod) {
    return (a >= b) ? (a - b) : (a + mod - b);
}
inline u32 mul_mod(u32 a, u32 b, u32 mod) {
    return (u64)a * b % mod;
}

inline u32 pow_mod(u32 a, u64 e, u32 mod) {
    u64 r = 1;
    u64 x = a;
    while (e) {
        if (e & 1) r = (r * x) % mod;
        x = (x * x) % mod;
        e >>= 1;
    }
    return (u32)r;
}
inline u32 inv_mod(u32 a, u32 mod) {
    // mod is prime
    return pow_mod(a, (u64)mod - 2, mod);
}

struct NTTCache {
    u32 mod = 0;
    u32 primitive_root = 0;
    int max_base = 0;

    u32 cache_n = 0;                 // power of two
    std::vector<u32> root;           // size cache_n: base^i
    std::vector<u32> iroot;          // size cache_n: inv(base)^i

    // inv_n_by_log[log] = inv(2^log) mod mod
    std::array<u32, 32> inv_n_by_log{};
    int inv_ready_upto = -1;         // max log filled
};

inline void build_root_table(NTTCache& C, u32 new_n) {
    // new_n must be power-of-two and <= 2^max_base
    if (new_n == 0 || (new_n & (new_n - 1)) != 0) {
        throw std::runtime_error("NTT cache: new_n must be power-of-two");
    }
    int logn = __builtin_ctz((unsigned)new_n);
    if (logn > C.max_base) {
        throw std::runtime_error("NTT cache: required N exceeds prime capability");
    }

    C.cache_n = new_n;
    C.root.assign(new_n, 0);
    C.iroot.assign(new_n, 0);

    // base = g^((mod-1)/N)
    const u32 base = pow_mod(C.primitive_root, (u64)(C.mod - 1) / (u64)new_n, C.mod);
    const u32 ibase = inv_mod(base, C.mod);

    C.root[0] = 1;
    for (u32 i = 1; i < new_n; ++i) {
        C.root[i] = mul_mod(C.root[i - 1], base, C.mod);
    }

    C.iroot[0] = 1;
    for (u32 i = 1; i < new_n; ++i) {
        C.iroot[i] = mul_mod(C.iroot[i - 1], ibase, C.mod);
    }

    // inv_n table (fill up to logn)
    for (int k = 0; k <= logn; ++k) {
        const u32 n_k = (u32)1u << k;
        C.inv_n_by_log[k] = inv_mod(n_k, C.mod);
    }
    C.inv_ready_upto = logn;
}

// Ensure cache has root tables for at least n (power of two)
inline NTTCache& ensure_ntt_cache(std::size_t prime_idx, u32 n) {
    static thread_local std::array<NTTCache, kNTT.size()> caches;

    NTTCache& C = caches[prime_idx];
    const auto prm = kNTT[prime_idx];

    if (C.mod == 0) {
        C.mod = prm.mod;
        C.primitive_root = prm.primitive_root;
        C.max_base = prm.max_base;
        C.cache_n = 0;
        C.inv_ready_upto = -1;
    }

    if (C.mod != prm.mod) {
        // should never happen
        throw std::runtime_error("NTT cache: modulus mismatch");
    }

    if (n == 0 || (n & (n - 1)) != 0) {
        throw std::runtime_error("NTT cache: n must be power-of-two");
    }

    // If cache_n is smaller than needed, rebuild to exactly n (or keep doubling if you prefer)
    if (C.cache_n < n) {
        build_root_table(C, n);
        return C;
    }

    // inv table might be too short if cache was built larger but inv_ready_upto not filled
    int logn = __builtin_ctz((unsigned)n);
    if (C.inv_ready_upto < logn) {
        for (int k = C.inv_ready_upto + 1; k <= logn; ++k) {
            const u32 n_k = (u32)1u << k;
            C.inv_n_by_log[k] = inv_mod(n_k, C.mod);
        }
        C.inv_ready_upto = logn;
    }

    return C;
}

// Cached NTT: no pow_mod in hot path
inline void ntt_cached(std::size_t prime_idx, std::vector<u32>& a, bool invert) {
    const u32 n = (u32)a.size();
    NTTCache& C = ensure_ntt_cache(prime_idx, n);
    const u32 mod = C.mod;

    // bit-reversal permutation (same as your old code, no extra memory)
    for (u32 i = 1, j = 0; i < n; ++i) {
        u32 bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // Use cached root table sized cache_n (>= n).
    // For stage len, wlen = g^((mod-1)/len) = base^(cache_n/len).
    const u32 cache_n = C.cache_n;

    for (u32 len = 2; len <= n; len <<= 1) {
        const u32 half = len >> 1;
        const u32 step = cache_n / len;
        const u32 wlen = invert ? C.iroot[step] : C.root[step];

        for (u32 i = 0; i < n; i += len) {
            u64 w = 1;
            for (u32 j = 0; j < half; ++j) {
                const u32 u = a[i + j];
                const u32 v = (u32)(w * a[i + j + half] % mod);
                a[i + j] = add_mod(u, v, mod);
                a[i + j + half] = sub_mod(u, v, mod);
                w = (w * wlen) % mod;
            }
        }
    }

    if (invert) {
        const int logn = __builtin_ctz((unsigned)n);
        const u32 inv_n = C.inv_n_by_log[logn];
        for (u32 i = 0; i < n; ++i) {
            a[i] = mul_mod(a[i], inv_n, mod);
        }
    }
}


// iterative Cooley–Tukey NTT (bit-reversal inside)
inline void ntt(std::vector<u32>& a, bool invert, u32 mod, u32 primitive_root) {
    const int n = (int)a.size();
    // bit reverse permutation
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        // wlen = g^((mod-1)/len)
        u32 wlen = pow_mod(primitive_root, (u64)(mod - 1) / (u64)len, mod);
        if (invert) wlen = inv_mod(wlen, mod);

        for (int i = 0; i < n; i += len) {
            u64 w = 1;
            const int half = len >> 1;
            for (int j = 0; j < half; ++j) {
                const u32 u = a[i + j];
                const u32 v = (u32)(w * a[i + j + half] % mod);
                a[i + j] = add_mod(u, v, mod);
                a[i + j + half] = sub_mod(u, v, mod);
                w = (w * wlen) % mod;
            }
        }
    }

    if (invert) {
        u32 inv_n = inv_mod((u32)n, mod);
        for (int i = 0; i < n; ++i) {
            a[i] = mul_mod(a[i], inv_n, mod);
        }
    }
}

inline std::vector<u32> convolution_mod_ntt(const std::vector<u32>& a,
                                            const std::vector<u32>& b,
                                            std::size_t prime_idx) {
    if (a.empty() || b.empty()) return {};

    const auto prm = kNTT[prime_idx];
    const u32 mod = prm.mod;

    std::size_t need = a.size() + b.size() - 1;
    std::size_t ntt_n = 1;
    while (ntt_n < need) ntt_n <<= 1;

    if ((int)__builtin_ctzll((unsigned long long)ntt_n) > prm.max_base) {
        throw std::runtime_error("convolution_mod_ntt: NTT length exceeds prime capability");
    }

    std::vector<u32> fa(ntt_n, 0), fb(ntt_n, 0);
    std::copy(a.begin(), a.end(), fa.begin());
    std::copy(b.begin(), b.end(), fb.begin());

    ntt_cached(prime_idx, fa, false);
    ntt_cached(prime_idx, fb, false);

    for (std::size_t i = 0; i < ntt_n; ++i) {
        fa[i] = (u32)((u64)fa[i] * fb[i] % mod);
    }

    ntt_cached(prime_idx, fa, true);

    fa.resize(need);
    return fa;
}


// Garner precomputation for CRT to target modulus p (u64)
// We will compute x mod p, assuming product(mods) > true integer convolution bound.
struct GarnerPrecomp {
    static constexpr int K = (int)kNTT.size();
    std::array<u32, K> mods{};
    std::array<u32, K> inv_coeff{}; // inv( prod_{t<i} mods[t] mod mods[i] ) mod mods[i]

    GarnerPrecomp() {
        for (int i = 0; i < K; ++i) mods[i] = kNTT[i].mod;

        for (int i = 0; i < K; ++i) {
            u32 mi = mods[i];
            u32 prod = 1;
            for (int t = 0; t < i; ++t) {
                prod = (u32)((u64)prod * (u64)(mods[t] % mi) % mi);
            }
            inv_coeff[i] = inv_mod(prod, mi);
        }
    }
};

inline const GarnerPrecomp& garner_precomp() {
    static const GarnerPrecomp pc;
    return pc;
}

// Compute CRT solution x modulo target prime p (u64), using Garner, without big integers.
// residues r[i] are mod mods[i] (u32).
inline u64 garner_crt_to_mod_p(const u32* r, u64 p) {
    const auto& pc = garner_precomp();
    constexpr int K = GarnerPrecomp::K;

    std::array<u64, K + 1> consts{};
    std::array<u64, K + 1> coeffs{};
    for (int i = 0; i <= K; ++i) {
        consts[i] = 0;
        coeffs[i] = 1;
    }

    for (int i = 0; i < K; ++i) {
        const u64 mi = pc.mods[i];

        u64 t = (u64)r[i];
        u64 cur = consts[i] % mi;
        // t = (r[i] - consts[i]) / coeffs[i]  (mod mi)
        t = (t + mi - cur) % mi;
        t = (t * (u64)pc.inv_coeff[i]) % mi;

        // update following moduli
        for (int j = i + 1; j < K; ++j) {
            const u64 mj = pc.mods[j];
            consts[j] = (consts[j] + (coeffs[j] % mj) * t) % mj;
            coeffs[j] = (coeffs[j] * mi) % mj;
        }
        // update target modulus p
        consts[K] = (consts[K] + (u128)(coeffs[K] % p) * (t % p) % p) % p;
        coeffs[K] = (u128)(coeffs[K] % p) * (mi % p) % p;
    }

    return consts[K] % p;
}

} // namespace detail


// 多项式：c_[i] 是 x^i 的系数，始终保持 trim 后（最高次系数非 0，除非零多项式）
class FpPoly {
public:
    using size_type = std::size_t;

    // 子乘积树：给 multipoint_eval / interpolation 预留
    struct SubproductTree;

private:
    // ---- fast div/mod helpers (internal) ----
    static FpPoly trunc_poly(const FpPoly& f, size_type k);
    static FpPoly reverse_poly(const FpPoly& f, size_type n);
    static FpPoly mul_trunc_poly(const FpPoly& a, const FpPoly& b, size_type k);
    static FpPoly inv_series_poly(const FpPoly& f, size_type k);

    FpPoly mod_slow_impl(const FpPoly& divisor) const;
    std::pair<FpPoly, FpPoly> divrem_slow_impl(const FpPoly& divisor) const;

    FpPoly mod_fast_impl(const FpPoly& divisor) const;
    std::pair<FpPoly, FpPoly> divrem_fast_impl(const FpPoly& divisor) const;

    const FpCtx* ctx_ = nullptr;
    std::vector<Fp> c_; // coefficients

    void require_ctx() const {
        if (!ctx_) throw std::logic_error("FpPoly: null ctx");
    }

    void require_compat(const FpPoly& other) const {
        require_ctx();
        other.require_ctx();
        if (ctx_->modulus() != other.ctx_->modulus()) {
            throw std::invalid_argument("FpPoly: modulus mismatch between polynomials");
        }
    }

public:
    // ---- constructors ----
    explicit FpPoly(const FpCtx& ctx) : ctx_(&ctx) {}

    FpPoly(const FpCtx& ctx, std::vector<Fp> coeffs)
        : ctx_(&ctx), c_(std::move(coeffs)) {
        // 确保每个系数都在 [0,p)（即使调用者传了非规范值也能工作）
        for (auto& x : c_) x.v %= ctx_->modulus();
        trim();
    }

    // u64 系数构造（自动 mod p）
    FpPoly(const FpCtx& ctx, std::initializer_list<u64> coeffs_u64)
        : ctx_(&ctx) {
        c_.reserve(coeffs_u64.size());
        for (u64 a : coeffs_u64) c_.push_back(ctx_->from_uint(a));
        trim();
    }

    // ---- basic access ----
    const FpCtx& ctx() const { require_ctx(); return *ctx_; }
    u64 modulus() const { return ctx().modulus(); }

    const std::vector<Fp>& coeffs() const noexcept { return c_; }
    std::vector<Fp>& coeffs_mut() noexcept { return c_; } // 谨慎使用：改完建议 trim()

    bool is_zero() const noexcept { return c_.empty(); }

    // degree: 零多项式返回 -1
    int degree() const noexcept {
        return c_.empty() ? -1 : static_cast<int>(c_.size()) - 1;
    }

    // 去掉最高位多余 0
    void trim() {
        while (!c_.empty() && c_.back().v == 0) {
            c_.pop_back();
        }
    }

    // 取第 i 项系数（越界视为 0）
    Fp coeff(size_type i) const noexcept {
        if (i >= c_.size()) return Fp{0};
        return c_[i];
    }

    // 常数项
    Fp constant_term() const noexcept {
        return c_.empty() ? Fp{0} : c_[0];
    }

    // 最高次系数（零多项式会抛异常）
    Fp leading_coeff() const {
        if (c_.empty()) throw std::domain_error("FpPoly::leading_coeff: zero polynomial");
        return c_.back();
    }

    // 设置系数：自动扩容 + 规约 + trim
    void set_coeff(size_type i, Fp value) {
        require_ctx();
        value.v %= ctx_->modulus();
        if (i >= c_.size()) c_.resize(i + 1, ctx_->zero());
        c_[i] = value;
        trim();
    }

    // ---- poly operations ----

    FpPoly add(const FpPoly& g) const {
        require_compat(g);
        const FpCtx& F = *ctx_;

        FpPoly r(F);
        const size_type n = std::max(c_.size(), g.c_.size());
        r.c_.assign(n, F.zero());

        for (size_type i = 0; i < n; ++i) {
            Fp a = (i < c_.size()) ? c_[i] : F.zero();
            Fp b = (i < g.c_.size()) ? g.c_[i] : F.zero();
            r.c_[i] = F.add(a, b);
        }
        r.trim();
        return r;
    }

    FpPoly sub(const FpPoly& g) const {
        require_compat(g);
        const FpCtx& F = *ctx_;

        FpPoly r(F);
        const size_type n = std::max(c_.size(), g.c_.size());
        r.c_.assign(n, F.zero());

        for (size_type i = 0; i < n; ++i) {
            Fp a = (i < c_.size()) ? c_[i] : F.zero();
            Fp b = (i < g.c_.size()) ? g.c_[i] : F.zero();
            r.c_[i] = F.sub(a, b);
        }
        r.trim();
        return r;
    }

    FpPoly mul(const FpPoly& g) const {
        require_compat(g);
        const FpCtx& F = *ctx_;

        if (is_zero() || g.is_zero()) return FpPoly(F);

        const std::size_t n = c_.size();
        const std::size_t m = g.c_.size();
        const std::size_t need = n + m - 1;

        // 小规模继续用朴素（避免 NTT 常数）
        // 你可以按机器再调阈值
        if (need <= 256 || std::min(n, m) <= 64) {
            // ---- naive O(n*m) ----
            std::vector<Fp> out(need, F.zero());
            for (std::size_t i = 0; i < n; ++i) {
                if (c_[i].v == 0) continue;
                for (std::size_t j = 0; j < m; ++j) {
                    if (g.c_[j].v == 0) continue;
                    out[i + j] = F.add(out[i + j], F.mul(c_[i], g.c_[j]));
                }
            }
            FpPoly r(F, std::move(out));
            r.trim();
            return r;
        }

        // ---- NTT (several friendly primes) + CRT/Garner back to modulus p ----
        // For each NTT prime, compute convolution residues, then Garner to mod p.
        const u64 p = F.modulus();

        // prepare inputs as u32 residues per prime
        std::vector<std::vector<detail::u32>> residues(detail::kNTT.size());

        // compute each modulus convolution
        for (std::size_t k = 0; k < detail::kNTT.size(); ++k) {
            const auto prm = detail::kNTT[k];

            std::vector<detail::u32> a(n), b(m);
            for (std::size_t i = 0; i < n; ++i) a[i] = (detail::u32)(c_[i].v % prm.mod);
            for (std::size_t j = 0; j < m; ++j) b[j] = (detail::u32)(g.c_[j].v % prm.mod);

            residues[k] = detail::convolution_mod_ntt(a, b, k);
            // residues[k].size() == need
        }

        // Garner combine to F_p
        std::vector<Fp> out(need, F.zero());
        std::vector<detail::u32> r(detail::kNTT.size());

        for (std::size_t idx = 0; idx < need; ++idx) {
            for (std::size_t k = 0; k < detail::kNTT.size(); ++k) {
                r[k] = residues[k][idx];
            }
            u64 val = detail::garner_crt_to_mod_p(r.data(), p);
            out[idx] = Fp{val};
        }

        FpPoly res(F, std::move(out));
        res.trim();
        return res;
    }


    // 标量乘：k * f(x)
    FpPoly scalar_mul(Fp k) const {
        require_ctx();
        const FpCtx& F = *ctx_;
        k.v %= F.modulus();

        if (k.v == 0 || is_zero()) return FpPoly(F);

        FpPoly r(F);
        r.c_.assign(c_.size(), F.zero());
        for (size_type i = 0; i < c_.size(); ++i) {
            r.c_[i] = F.mul(c_[i], k);
        }
        r.trim();
        return r;
    }

    // 导数
    FpPoly derivative() const {
        require_ctx();
        const FpCtx& F = *ctx_;

        if (c_.size() <= 1) return FpPoly(F);

        FpPoly r(F);
        r.c_.assign(c_.size() - 1, F.zero());

        for (size_type i = 1; i < c_.size(); ++i) {
            // (a_i * i) x^{i-1}
            Fp ii = F.from_uint(static_cast<u64>(i)); // 自动 i mod p（特征 p 的情况也正确）
            r.c_[i - 1] = F.mul(c_[i], ii);
        }
        r.trim();
        return r;
    }

    // 单点求值（Horner）
    Fp eval(Fp x) const {
        require_ctx();
        const FpCtx& F = *ctx_;
        x.v %= F.modulus();

        Fp acc = F.zero();
        for (size_type i = c_.size(); i-- > 0;) {
            acc = F.add(F.mul(acc, x), c_[i]);
        }
        return acc;
    }

    // ---- division/remainder: 为 multipoint_eval（树）提供基础 ----
    // Long division: this = q*div + r
    std::pair<FpPoly, FpPoly> divrem(const FpPoly& divisor) const {
        require_compat(divisor);
        const FpCtx& F = *ctx_;

        if (divisor.is_zero()) {
            throw std::domain_error("FpPoly::divrem: division by zero polynomial");
        }
        if (this->is_zero()) {
            return {FpPoly(F), FpPoly(F)};
        }

        const int degA = this->degree();
        const int degB = divisor.degree();
        if (degA < degB) {
            return {FpPoly(F), *this};
        }
        if (degB == 0) {
            Fp invb = F.inv(divisor.c_[0]);
            return {this->scalar_mul(invb), FpPoly(F)};
        }

        const size_type n = static_cast<size_type>(degA + 1);
        const size_type m = static_cast<size_type>(degB + 1);

        // 小规模用 slow，常数更低
        if (n <= 512 || m <= 64 || (n - m) <= 64) {
            return divrem_slow_impl(divisor);
        }
        return divrem_fast_impl(divisor);
    }


    FpPoly mod(const FpPoly& divisor) const {
        require_compat(divisor);
        const FpCtx& F = *ctx_;

        if (divisor.is_zero()) {
            throw std::domain_error("FpPoly::mod: division by zero polynomial");
        }
        if (this->is_zero()) return FpPoly(F);

        const int degA = this->degree();
        const int degB = divisor.degree();
        if (degA < degB) return *this;
        if (degB == 0) return FpPoly(F);

        const size_type n = static_cast<size_type>(degA + 1);
        const size_type m = static_cast<size_type>(degB + 1);

        if (n <= 512 || m <= 64 || (n - m) <= 64) {
            return mod_slow_impl(divisor);
        }
        return mod_fast_impl(divisor);
    }



    // ---- convenience operations & operators ----

    // 系数逐项取负
    FpPoly neg_poly() const {
        require_ctx();
        const FpCtx& F = *ctx_;
        if (is_zero()) return FpPoly(F);

        FpPoly r(F);
        r.c_.assign(c_.size(), F.zero());
        for (size_type i = 0; i < c_.size(); ++i) {
            r.c_[i] = F.neg(c_[i]);
        }
        r.trim();
        return r;
    }

    // 忽略尾部 0 的“逻辑相等”（不修改对象）
    bool equals(const FpPoly& g) const noexcept {
        if (!ctx_ || !g.ctx_) return false;
        if (ctx_->modulus() != g.ctx_->modulus()) return false;

        size_type na = c_.size();
        while (na > 0 && c_[na - 1].v == 0) --na;

        size_type nb = g.c_.size();
        while (nb > 0 && g.c_[nb - 1].v == 0) --nb;

        if (na != nb) return false;
        for (size_type i = 0; i < na; ++i) {
            if (c_[i].v != g.c_[i].v) return false;
        }
        return true;
    }

    // 运算符：正确优先（内部直接调用已有成员函数）
    FpPoly operator+(const FpPoly& g) const { return add(g); }
    FpPoly operator-(const FpPoly& g) const { return sub(g); }
    FpPoly operator*(const FpPoly& g) const { return mul(g); }
    FpPoly operator-() const { return neg_poly(); }

    FpPoly& operator+=(const FpPoly& g) { *this = add(g); return *this; }
    FpPoly& operator-=(const FpPoly& g) { *this = sub(g); return *this; }
    FpPoly& operator*=(const FpPoly& g) { *this = mul(g); return *this; }

    // 标量乘
    FpPoly operator*(Fp k) const { return scalar_mul(k); }
    FpPoly& operator*=(Fp k) { *this = scalar_mul(k); return *this; }
    friend FpPoly operator*(Fp k, const FpPoly& f) { return f.scalar_mul(k); }

    // 求值便捷写法：f(x)
    Fp operator()(Fp x) const { return eval(x); }

    // 逻辑相等
    bool operator==(const FpPoly& g) const noexcept { return equals(g); }
    bool operator!=(const FpPoly& g) const noexcept { return !equals(g); }

    // ---- remainder tree (reserved + useful) ----
    // rem[level][idx] = f mod tree.levels[level][idx]
    // level=0 是叶子 (x-x_i)，此时 rem[0][i] 的常数项就是 f(x_i)
    std::vector<std::vector<FpPoly>> remainder_tree(const SubproductTree& tree) const;


    // ---- multipoint evaluation ----

    // 朴素多点求值：O(n*deg)
    std::vector<Fp> multipoint_eval_naive(const std::vector<Fp>& xs) const {
        require_ctx();
        const FpCtx& F = *ctx_;
        std::vector<Fp> ys;
        ys.reserve(xs.size());
        for (Fp x : xs) {
            x.v %= F.modulus();
            ys.push_back(eval(x));
        }
        return ys;
    }

    // 用子乘积树做多点求值（正确，但当前不追求性能）
    std::vector<Fp> multipoint_eval(const SubproductTree& tree) const;

    // 便捷：内部建树再求值（你也可以后续加策略：小 n 用 naive，大 n 用 tree）
    std::vector<Fp> multipoint_eval(const std::vector<Fp>& xs) const;

    // ---- interpolation interfaces (reserved) ----

    // O(n^2) Lagrange 插值：用于基准正确性（可直接用）
    static FpPoly interpolate_lagrange_naive(const FpCtx& ctx,
                                            const std::vector<Fp>& xs,
                                            const std::vector<Fp>& ys);

    // 预留：基于子乘积树的插值（TODO）
    static FpPoly interpolate_subproduct_tree(const SubproductTree& tree,
                                             const std::vector<Fp>& ys);
};

inline pf::FpPoly pf::FpPoly::trunc_poly(const FpPoly& f, size_type k) {
    f.require_ctx();
    const FpCtx& F = *f.ctx_;
    if (k == 0 || f.c_.empty()) return FpPoly(F);
    const size_type take = std::min(k, f.c_.size());
    std::vector<Fp> v(take);
    for (size_type i = 0; i < take; ++i) v[i] = f.c_[i];
    return FpPoly(F, std::move(v));
}

inline pf::FpPoly pf::FpPoly::reverse_poly(const FpPoly& f, size_type n) {
    f.require_ctx();
    const FpCtx& F = *f.ctx_;
    if (n == 0) return FpPoly(F);
    std::vector<Fp> v(n, F.zero());
    for (size_type i = 0; i < n; ++i) {
        const size_type src = n - 1 - i;
        if (src < f.c_.size()) v[i] = f.c_[src];
    }
    return FpPoly(F, std::move(v));
}

inline pf::FpPoly pf::FpPoly::mul_trunc_poly(const FpPoly& a, const FpPoly& b, size_type k) {
    a.require_compat(b);
    const FpCtx& F = *a.ctx_;
    if (k == 0) return FpPoly(F);
    if (a.is_zero() || b.is_zero()) return FpPoly(F);

    FpPoly prod = a.mul(b);
    if (prod.c_.size() > k) prod.c_.resize(k);
    prod.trim();
    return prod;
}

inline pf::FpPoly pf::FpPoly::inv_series_poly(const FpPoly& f, size_type k) {
    f.require_ctx();
    const FpCtx& F = *f.ctx_;
    if (k == 0) return FpPoly(F);
    if (f.c_.empty() || f.c_[0].v == 0) {
        throw std::domain_error("inv_series_poly: f[0] must be non-zero");
    }

    FpPoly g(F, { F.inv(f.c_[0]) });

    size_type cur = 1;
    const Fp two = F.from_uint(2);

    while (cur < k) {
        const size_type nxt = std::min(cur * 2, k);

        FpPoly f_tr = trunc_poly(f, nxt);
        FpPoly t = mul_trunc_poly(f_tr, g, nxt);

        std::vector<Fp> uc(nxt, F.zero());
        for (size_type i = 0; i < nxt; ++i) {
            Fp ti = (i < t.c_.size()) ? t.c_[i] : F.zero();
            uc[i] = F.neg(ti);
        }
        uc[0] = F.add(uc[0], two);
        FpPoly u(F, std::move(uc));

        g = mul_trunc_poly(g, u, nxt);
        cur = nxt;
    }

    return trunc_poly(g, k);
}

inline std::pair<pf::FpPoly, pf::FpPoly>
pf::FpPoly::divrem_slow_impl(const FpPoly& divisor) const {
    require_compat(divisor);
    const FpCtx& F = *ctx_;

    if (divisor.is_zero()) throw std::domain_error("divrem_slow: div by zero");
    if (this->is_zero()) return {FpPoly(F), FpPoly(F)};

    const int degA = this->degree();
    const int degB = divisor.degree();
    if (degA < degB) return {FpPoly(F), *this};

    if (degB == 0) {
        Fp invb = F.inv(divisor.c_[0]);
        return {this->scalar_mul(invb), FpPoly(F)};
    }

    std::vector<Fp> r = c_;
    std::vector<Fp> q(static_cast<size_type>(degA - degB + 1), F.zero());

    const Fp lcB = divisor.c_.back();
    const bool monic = (lcB.v == 1);
    const Fp inv_lcB = monic ? F.one() : F.inv(lcB);

    int degR = degA;
    while (degR >= degB) {
        const Fp lead = r[static_cast<size_type>(degR)];
        if (lead.v != 0) {
            const Fp factor = monic ? lead : F.mul(lead, inv_lcB);
            const int k = degR - degB;
            q[static_cast<size_type>(k)] = factor;

            for (int i = 0; i < degB; ++i) {
                const Fp di = divisor.c_[static_cast<size_type>(i)];
                if (di.v == 0) continue;
                const size_type idx = static_cast<size_type>(i + k);
                r[idx] = F.sub(r[idx], F.mul(factor, di));
            }
            r[static_cast<size_type>(degR)] = F.zero();
        } else {
            r[static_cast<size_type>(degR)] = F.zero();
        }

        --degR;
        while (degR >= 0 && r[static_cast<size_type>(degR)].v == 0) --degR;
    }

    FpPoly Q(F, std::move(q));
    if (degR < 0) return {Q, FpPoly(F)};

    r.resize(static_cast<size_type>(degR + 1));
    FpPoly R(F, std::move(r));
    return {Q, R};
}

inline pf::FpPoly pf::FpPoly::mod_slow_impl(const FpPoly& divisor) const {
    return divrem_slow_impl(divisor).second;
}

inline std::pair<pf::FpPoly, pf::FpPoly>
pf::FpPoly::divrem_fast_impl(const FpPoly& divisor) const {
    require_compat(divisor);
    const FpCtx& F = *ctx_;

    if (divisor.is_zero()) throw std::domain_error("divrem_fast: div by zero");
    if (this->is_zero()) return {FpPoly(F), FpPoly(F)};

    const int degA = this->degree();
    const int degB = divisor.degree();
    if (degA < degB) return {FpPoly(F), *this};

    if (degB == 0) {
        Fp invb = F.inv(divisor.c_[0]);
        return {this->scalar_mul(invb), FpPoly(F)};
    }

    const size_type n = static_cast<size_type>(degA + 1);
    const size_type m = static_cast<size_type>(degB + 1);
    const size_type k = n - m + 1;

    FpPoly Ar = reverse_poly(*this, n);
    FpPoly Br = reverse_poly(divisor, m);

    FpPoly Br_inv = inv_series_poly(Br, k);

    FpPoly Ar_tr = trunc_poly(Ar, k);
    FpPoly qrev = mul_trunc_poly(Ar_tr, Br_inv, k);

    FpPoly q = reverse_poly(qrev, k);
    q.trim();

    FpPoly r = this->sub(divisor.mul(q));
    r.trim();

    return {q, r};
}

inline pf::FpPoly pf::FpPoly::mod_fast_impl(const FpPoly& divisor) const {
    return divrem_fast_impl(divisor).second;
}


// 子乘积树：levels[0] = (x - x_i)，levels.back()[0] = Π(x - x_i)
struct FpPoly::SubproductTree {
    const FpCtx* ctx = nullptr;
    std::vector<Fp> points;
    std::vector<std::vector<FpPoly>> levels;

    SubproductTree() = default;
    explicit SubproductTree(const FpCtx& c) : ctx(&c) {}

    bool empty() const noexcept { return points.empty(); }
    std::size_t n_points() const noexcept { return points.size(); }
    std::size_t n_levels() const noexcept { return levels.size(); }

    const FpPoly& root() const {
        if (levels.empty() || levels.back().empty()) {
            throw std::logic_error("SubproductTree::root: empty tree");
        }
        return levels.back()[0];
    }

    // 构建子乘积树（正确优先，乘法用朴素 mul）
    static SubproductTree build(const FpCtx& ctx, const std::vector<Fp>& xs) {
        SubproductTree T(ctx);
        T.points = xs;
        for (auto& x : T.points) x.v %= ctx.modulus();

        if (T.points.empty()) return T;

        // level 0
        std::vector<FpPoly> level0;
        level0.reserve(T.points.size());
        for (const auto& xi : T.points) {
            // (x - xi) = (-xi) + 1*x
            FpPoly leaf(ctx);
            leaf.coeffs_mut().reserve(2);
            leaf.coeffs_mut().push_back(ctx.neg(xi));
            leaf.coeffs_mut().push_back(ctx.one());
            leaf.trim();
            level0.push_back(std::move(leaf));
        }
        T.levels.push_back(std::move(level0));

        // upper levels
        while (T.levels.back().size() > 1) {
            const auto& prev = T.levels.back();
            std::vector<FpPoly> nxt;
            nxt.reserve((prev.size() + 1) / 2);

            for (std::size_t i = 0; i < prev.size(); i += 2) {
                if (i + 1 < prev.size()) {
                    nxt.push_back(prev[i].mul(prev[i + 1]));
                } else {
                    nxt.push_back(prev[i]); // carry
                }
            }
            T.levels.push_back(std::move(nxt));
        }

        return T;
    }
};

inline std::vector<std::vector<FpPoly>>
FpPoly::remainder_tree(const SubproductTree& tree) const {
    require_ctx();
    if (!tree.ctx) {
        throw std::invalid_argument("FpPoly::remainder_tree: tree.ctx is null");
    }
    if (ctx_->modulus() != tree.ctx->modulus()) {
        throw std::invalid_argument("FpPoly::remainder_tree: modulus mismatch");
    }

    const FpCtx& F = *ctx_;

    const std::size_t n = tree.n_points();
    if (n == 0) return {}; // 空树

    if (tree.levels.empty() || tree.levels[0].size() != n) {
        throw std::invalid_argument("FpPoly::remainder_tree: malformed tree (levels[0] size mismatch)");
    }

    const std::size_t L = tree.n_levels();
    if (L == 0) {
        throw std::invalid_argument("FpPoly::remainder_tree: malformed tree (no levels)");
    }

    // rem 与 tree.levels 形状一致
    std::vector<std::vector<FpPoly>> rem;
    rem.reserve(L);
    for (std::size_t level = 0; level < L; ++level) {
        rem.emplace_back(tree.levels[level].size(), FpPoly(F));
    }

    const std::size_t top = L - 1;

    // 顶层：f mod root
    rem[top][0] = this->mod(tree.levels[top][0]);

    // 自顶向下传播：child_rem = parent_rem mod child_modulus_poly
    for (std::size_t level = top; level-- > 0;) {
        // 从 level+1 的每个结点，传播到 level 的孩子
        const std::size_t parent_level = level + 1;

        for (std::size_t idx = 0; idx < rem[parent_level].size(); ++idx) {
            const FpPoly& r_parent = rem[parent_level][idx];

            const std::size_t left = idx * 2;
            if (left < rem[level].size()) {
                rem[level][left] = r_parent.mod(tree.levels[level][left]);
            }

            const std::size_t right = left + 1;
            if (right < rem[level].size()) {
                rem[level][right] = r_parent.mod(tree.levels[level][right]);
            }
        }
    }

    return rem;
}


// ---- multipoint_eval implementation ----
inline std::vector<Fp> FpPoly::multipoint_eval(const SubproductTree& tree) const {
    require_ctx();
    if (!tree.ctx) {
        throw std::invalid_argument("FpPoly::multipoint_eval(tree): tree.ctx is null");
    }
    if (ctx_->modulus() != tree.ctx->modulus()) {
        throw std::invalid_argument("FpPoly::multipoint_eval(tree): modulus mismatch");
    }

    const FpCtx& F = *ctx_;
    const std::size_t n = tree.n_points();
    std::vector<Fp> ys(n, F.zero());
    if (n == 0) return ys;

    auto rem = remainder_tree(tree);

    // rem[0][i] = f mod (x-x_i) = 常数多项式 [f(x_i)]
    if (rem.empty() || rem[0].size() != n) {
        throw std::logic_error("FpPoly::multipoint_eval: unexpected remainder_tree shape");
    }

    for (std::size_t i = 0; i < n; ++i) {
        ys[i] = rem[0][i].constant_term();
    }
    return ys;
}


inline std::vector<Fp> FpPoly::multipoint_eval(const std::vector<Fp>& xs) const {
    require_ctx();
    SubproductTree T = SubproductTree::build(*ctx_, xs);
    return multipoint_eval(T);
}

// ---- interpolation ----
inline FpPoly FpPoly::interpolate_lagrange_naive(const FpCtx& ctx,
                                                 const std::vector<Fp>& xs,
                                                 const std::vector<Fp>& ys) {
    if (xs.size() != ys.size()) {
        throw std::invalid_argument("interpolate_lagrange_naive: xs.size != ys.size");
    }
    const std::size_t n = xs.size();
    FpPoly result(ctx);
    if (n == 0) return result;

    // 归一化输入
    std::vector<Fp> X = xs, Y = ys;
    for (auto& x : X) x.v %= ctx.modulus();
    for (auto& y : Y) y.v %= ctx.modulus();

    // result = Σ y_i * Π_{j!=i} (x-x_j)/(x_i-x_j)
    for (std::size_t i = 0; i < n; ++i) {
        FpPoly numer(ctx, {1}); // 1
        Fp denom = ctx.one();

        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;

            // numer *= (x - x_j)
            FpPoly lin(ctx);
            lin.coeffs_mut().reserve(2);
            lin.coeffs_mut().push_back(ctx.neg(X[j]));
            lin.coeffs_mut().push_back(ctx.one());
            lin.trim();
            numer = numer.mul(lin);

            // denom *= (x_i - x_j)
            denom = ctx.mul(denom, ctx.sub(X[i], X[j]));
        }

        Fp scale = ctx.div(Y[i], denom); // y_i / denom
        result = result.add(numer.scalar_mul(scale));
    }

    result.trim();
    return result;
}

inline FpPoly FpPoly::interpolate_subproduct_tree(const SubproductTree& tree,
                                                  const std::vector<Fp>& ys) {
    if (!tree.ctx) {
        throw std::invalid_argument("interpolate_subproduct_tree: tree.ctx is null");
    }
    const FpCtx& F = *tree.ctx;

    const std::size_t n = tree.n_points();
    if (ys.size() != n) {
        throw std::invalid_argument("interpolate_subproduct_tree: ys.size != number of points");
    }
    if (n == 0) {
        return FpPoly(F); // 零多项式
    }
    if (tree.levels.empty() || tree.levels[0].size() != n) {
        throw std::invalid_argument("interpolate_subproduct_tree: malformed tree (levels[0] size mismatch)");
    }

    // 1) G(x) = Π (x - x_i)
    const FpPoly& G = tree.root();

    // 2) dG = G'(x)
    FpPoly dG = G.derivative();

    // 3) 计算 dG(x_i)
    std::vector<Fp> dvals = dG.multipoint_eval(tree);
    if (dvals.size() != n) {
        throw std::logic_error("interpolate_subproduct_tree: unexpected dvals size");
    }

    // 4) 一次性求逆：dvals[i] <- inv(dG(x_i))
    for (std::size_t i = 0; i < n; ++i) {
        if (dvals[i].v == 0) {
            // 通常意味着 x_i 不互异导致 dG(x_i)=0
            throw std::domain_error("interpolate_subproduct_tree: dG(x_i)=0 (points likely not distinct)");
        }
    }
    F.batch_inv(dvals);

    // 5) a_i = y_i * inv(dG(x_i))
    std::vector<Fp> a(n, F.zero());
    for (std::size_t i = 0; i < n; ++i) {
        Fp yi = ys[i];
        yi.v %= F.modulus();
        a[i] = F.mul(yi, dvals[i]);
    }


    // 6) 自底向上合并：
    // 在每个结点 S 上维护 F_S(x) = Σ_{i in S} a_i * (M_S(x)/(x-x_i))
    // 叶子：F_{ {i} } = a_i (常数多项式)
    std::vector<std::vector<FpPoly>> Flevels;
    Flevels.reserve(tree.n_levels());

    // level 0 (leaf): constant polys
    {
        std::vector<FpPoly> L0;
        L0.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            FpPoly leaf(F);
            if (a[i].v != 0) {
                leaf.coeffs_mut().push_back(a[i]); // 常数项
            }
            leaf.trim();
            L0.push_back(std::move(leaf));
        }
        Flevels.push_back(std::move(L0));
    }

    // level k>0: combine pairs
    for (std::size_t level = 1; level < tree.levels.size(); ++level) {
        const auto& prevF = Flevels[level - 1];
        const auto& prevM = tree.levels[level - 1]; // 同层的 M 子树多项式（产品）

        if (prevF.size() != prevM.size()) {
            throw std::invalid_argument("interpolate_subproduct_tree: tree malformed (prevF size != prevM size)");
        }

        std::vector<FpPoly> cur;
        cur.reserve((prevF.size() + 1) / 2);

        for (std::size_t i = 0; i < prevF.size(); i += 2) {
            if (i + 1 < prevF.size()) {
                // parent = F_left*M_right + F_right*M_left
                FpPoly t1 = prevF[i].mul(prevM[i + 1]);
                FpPoly t2 = prevF[i + 1].mul(prevM[i]);
                FpPoly sum = t1.add(t2);
                sum.trim();
                cur.push_back(std::move(sum));
            } else {
                // odd carry: parent is same as child
                cur.push_back(prevF[i]);
            }
        }

        Flevels.push_back(std::move(cur));
    }

    if (Flevels.empty() || Flevels.back().empty()) {
        throw std::logic_error("interpolate_subproduct_tree: internal error (empty result)");
    }

    FpPoly result = Flevels.back()[0];
    result.trim();
    return result;
}


// 简单打印（调试用）
inline std::ostream& operator<<(std::ostream& os, const FpPoly& f) {
    if (f.is_zero()) return os << "0";
    bool first = true;
    for (std::size_t i = 0; i < f.coeffs().size(); ++i) {
        Fp ci = f.coeffs()[i];
        if (ci.v == 0) continue;
        if (!first) os << " + ";
        first = false;
        os << ci.v;
        if (i >= 1) os << "*x";
        if (i >= 2) os << "^" << i;
    }
    return os;
}

} // namespace pf

#endif // FP_POLY_HPP
