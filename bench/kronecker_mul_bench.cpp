#include <poly_interp/prime_field.hpp>
#include <poly_interp/fp_poly.hpp>

#include <gmp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using pf::Fp;
using pf::FpCtx;
using pf::FpPoly;
using pf::u64;

// ------------------------------ helpers ------------------------------
static u64 mersenne31() { return (1ULL << 31) - 1; }
static u64 mersenne61() { return (1ULL << 61) - 1; }

static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && std::memcmp(s.data(), prefix.data(), prefix.size()) == 0;
}

static u64 parse_u64(const std::string& s) {
    std::size_t pos = 0;
    unsigned long long v = std::stoull(s, &pos, 10);
    if (pos != s.size()) throw std::invalid_argument("invalid integer: " + s);
    return static_cast<u64>(v);
}

static u64 mix64(u64 h, u64 x) {
    // Simple FNV-1a style mix to keep results live
    h ^= x;
    h *= 1099511628211ULL;
    h ^= h >> 32;
    return h;
}

static u64 poly_fingerprint(const FpPoly& f) {
    u64 h = 1469598103934665603ULL;
    int deg = f.degree();
    h = mix64(h, static_cast<u64>(deg));
    if (deg < 0) return h;

    auto pick = [&](std::size_t i) {
        Fp ci = f.coeff(i);
        h = mix64(h, ci.v);
    };

    pick(0);
    if (deg >= 1) pick(1);
    if (deg >= 2) pick(2);
    std::size_t mid = static_cast<std::size_t>(deg) / 2;
    pick(mid);
    h = mix64(h, f.leading_coeff().v);
    return h;
}

using u128 = unsigned __int128;

// ------------------------------ Kronecker mul ------------------------------
static unsigned choose_base_bits(u64 modulus, std::size_t min_len) {
    // Upper bound: max coeff <= min_len * (p-1)^2. Take a couple safety bits.
    double lp = std::log2(static_cast<double>(modulus));
    double bound = 2.0 * lp + std::log2(static_cast<double>(std::max<std::size_t>(min_len, 1))) + 2.0;
    unsigned bits = static_cast<unsigned>(std::ceil(bound));
    if (bits < 16) bits = 16;
    return bits;
}

// mask helper that avoids UB on shift==64
static inline mp_limb_t mask_bits(unsigned bits) {
    if (bits == 0) return 0;
    if (bits >= GMP_NUMB_BITS) return ~(mp_limb_t)0;
    return ((mp_limb_t)1 << bits) - 1;
}

// Pack coefficients into contiguous bits (base 2^base_bits) and import once.
// Thread-local scratch to reuse buffers and pow2 tables during pack/unpack.
struct PackScratch {
    std::vector<u64> pow2; // 2^k mod p for current base_bits/mod
    unsigned pow_bits = 0;
    u64 pow_mod = 0;

    void ensure_pow2(unsigned bits, u64 mod) {
        if (pow_bits == bits && pow_mod == mod && !pow2.empty()) return;
        pow_bits = bits;
        pow_mod = mod;
        pow2.assign(bits + 1, 0);
        pow2[0] = 1 % mod;
        for (unsigned i = 1; i <= bits; ++i) {
            pow2[i] = (u64)((u128)pow2[i - 1] * 2 % mod);
        }
    }
};

struct MpnScratch {
    std::vector<mp_limb_t> A;
    std::vector<mp_limb_t> B;
    std::vector<mp_limb_t> C;
};

static mp_size_t pack_poly_bits(std::vector<mp_limb_t>& out,
                                const std::vector<Fp>& coeffs,
                                unsigned base_bits) {
    if (coeffs.empty()) {
        out.clear();
        return 0;
    }

    const unsigned limb_bits = GMP_NUMB_BITS;
    const std::size_t total_bits = (std::size_t)coeffs.size() * (std::size_t)base_bits;
    const std::size_t limb_count = (total_bits + limb_bits - 1) / limb_bits;

    out.assign(limb_count, 0);

    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        u128 val = (u128)coeffs[i].v;
        std::size_t bitpos = (std::size_t)i * (std::size_t)base_bits;
        std::size_t idx = bitpos / limb_bits;
        unsigned offset = (unsigned)(bitpos % limb_bits);

        unsigned consumed = 0;
        unsigned bits_left = base_bits;
        while (bits_left > 0) {
            const unsigned take = std::min<unsigned>(bits_left, limb_bits - offset);
            mp_limb_t chunk = 0;
            if (consumed < 128) { // u128 width guard
                const unsigned safe_take = std::min<unsigned>(take, 128 - consumed);
                chunk = (mp_limb_t)((val >> consumed) & mask_bits(safe_take));
            }
            out[idx] |= (chunk << offset);

            consumed += take;
            bits_left -= take;
            ++idx;
            offset = 0;
        }
    }

    mp_size_t used = (mp_size_t)limb_count;
    while (used > 0 && out[(std::size_t)used - 1] == 0) --used;
    if (used == 0) used = 1; // keep at least one limb live
    return used;
}

// Export bits and recover coefficients modulo p.
static void unpack_poly_bits(std::vector<Fp>& out,
                             const mp_limb_t* limbs,
                             mp_size_t limb_count,
                             unsigned base_bits,
                             const FpCtx& F,
                             PackScratch& S) {
    const unsigned limb_bits = GMP_NUMB_BITS;
    const std::size_t limb_count_s = (limb_count > 0) ? (std::size_t)limb_count : 0;

    if (limb_count_s == 0 || limbs == nullptr) {
        for (auto& x : out) x = F.zero();
        return;
    }

    // Precompute 2^k mod p up to base_bits.
    const u64 mod = F.modulus();
    S.ensure_pow2(base_bits, mod);

    for (std::size_t i = 0; i < out.size(); ++i) {
        const std::size_t bitpos = (std::size_t)i * (std::size_t)base_bits;
        std::size_t idx = bitpos / limb_bits;
        if (idx >= limb_count_s) {
            out[i] = F.zero();
            continue;
        }
        unsigned offset = (unsigned)(bitpos % limb_bits);

        unsigned collected = 0;
        u64 acc_mod = 0;
        while (collected < base_bits && idx < limb_count_s) {
            const unsigned take = std::min<unsigned>(base_bits - collected, limb_bits - offset);
            const mp_limb_t limb = limbs[idx];
            mp_limb_t chunk = limb >> offset;
            if (take < GMP_NUMB_BITS) {
                chunk &= mask_bits(take);
            }

            // acc_mod = acc_mod + chunk * 2^{collected} (mod p)
            const u64 term = (u64)((u128)(chunk % mod) * (u128)S.pow2[collected] % mod);
            acc_mod += term;
            if (acc_mod >= mod) acc_mod -= mod;

            collected += take;
            ++idx;
            offset = 0;
        }

        out[i] = F.from_uint(acc_mod);
    }
}

static FpPoly mul_kronecker(const FpPoly& a, const FpPoly& b, unsigned* used_base_bits = nullptr) {
    const FpCtx& F = a.ctx();
    if (a.modulus() != b.modulus()) {
        throw std::invalid_argument("mul_kronecker: modulus mismatch");
    }

    if (a.is_zero() || b.is_zero()) return FpPoly(F);

    const std::size_t n = a.coeffs().size();
    const std::size_t m = b.coeffs().size();
    const std::size_t min_len = std::min(n, m);
    const unsigned base_bits = choose_base_bits(F.modulus(), min_len);
    if (used_base_bits) *used_base_bits = base_bits;

    static thread_local MpnScratch S;
    static thread_local PackScratch PS;

    const mp_size_t limbs_a = pack_poly_bits(S.A, a.coeffs(), base_bits);
    const mp_size_t limbs_b = pack_poly_bits(S.B, b.coeffs(), base_bits);

    const mp_size_t limbs_res = limbs_a + limbs_b;
    S.C.resize((std::size_t)limbs_res);

    if (limbs_a >= limbs_b) {
        mpn_mul(S.C.data(), S.A.data(), limbs_a, S.B.data(), limbs_b);
    } else {
        mpn_mul(S.C.data(), S.B.data(), limbs_b, S.A.data(), limbs_a);
    }

    mp_size_t used_limbs = limbs_res;
    while (used_limbs > 0 && S.C[(std::size_t)used_limbs - 1] == 0) --used_limbs;
    if (used_limbs == 0) used_limbs = 1;

    std::vector<Fp> out(n + m - 1, F.zero());
    unpack_poly_bits(out, S.C.data(), used_limbs, base_bits, F, PS);

    FpPoly res(F, std::move(out));
    res.trim();
    return res;
}

// ------------------------------ bench harness ------------------------------
struct Timings {
    double avg_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
    u64 fingerprint = 0;
};

static FpPoly random_poly(const FpCtx& F,
                          std::size_t len,
                          std::mt19937_64& rng,
                          std::uniform_int_distribution<u64>& dist_nonzero) {
    std::vector<Fp> c(len);
    for (std::size_t i = 0; i < len; ++i) c[i] = Fp{dist_nonzero(rng)};
    // distribution already excludes zero; guard in case of future changes
    if (len > 0 && c.back().v == 0) c.back().v = 1;
    FpPoly p(F, std::move(c));
    return p;
}

static void warmup_baseline_mul(const FpPoly& a, const FpPoly& b, int iters) {
    volatile u64 sink = 0;
    for (int i = 0; i < iters; ++i) {
        FpPoly r = a.mul(b);
        sink ^= poly_fingerprint(r);
    }
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
}

static void warmup_kron_mul(const FpPoly& a, const FpPoly& b, int iters) {
    volatile u64 sink = 0;
    for (int i = 0; i < iters; ++i) {
        FpPoly r = mul_kronecker(a, b);
        sink ^= poly_fingerprint(r);
    }
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
}

static Timings bench_mul_baseline(const std::vector<FpPoly>& a,
                                  const std::vector<FpPoly>& b) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(a.size());

    if (repeats > 0) warmup_baseline_mul(a[0], b[0], 1);

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        FpPoly r = a[rep].mul(b[rep]);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(r), (u64)rep + 1);
        sink ^= t.fingerprint;
    }
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_mul_kronecker(const std::vector<FpPoly>& a,
                                   const std::vector<FpPoly>& b,
                                   unsigned& base_bits_out) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(a.size());

    if (repeats > 0) warmup_kron_mul(a[0], b[0], 1);

    unsigned bb = 0;
    for (int rep = 0; rep < repeats; ++rep) {
        unsigned used_bits = 0;
        auto t0 = std::chrono::steady_clock::now();
        FpPoly r = mul_kronecker(a[rep], b[rep], &used_bits);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(r), (u64)rep + 1);
        sink ^= t.fingerprint;
        bb = std::max(bb, used_bits);
    }
    base_bits_out = bb;
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

struct Options {
    std::string prime_mode = "61"; // 31 | 61
    int min_pow = 12;
    int max_pow = 16;
    int repeats = 3;
    u64 seed = 1234567;
    bool csv = false;
    bool verify = true;
};

static void print_help(const char* argv0) {
    std::cout
        << "Kronecker multiplication experiment\n"
        << "Usage:\n"
        << "  " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --prime=31|61        choose Mersenne prime (default: 61)\n"
        << "  --min_pow=<k>        min k where n=2^k (default: 12)\n"
        << "  --max_pow=<k>        max k where n=2^k (default: 16)\n"
        << "  --repeats=<r>        repeats per n (default: 2)\n"
        << "  --seed=<s>           RNG seed (default: 1234567)\n"
        << "  --csv=0|1            CSV output (default: 0)\n"
        << "  --verify=0|1         check kronecker result against baseline (default: 1)\n"
        << "  --help               show this help\n";
}

static Options parse_args(int argc, char** argv) {
    Options opt;

    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        if (a == "--help") {
            print_help(argv[0]);
            std::exit(0);
        } else if (starts_with(a, "--prime=")) {
            opt.prime_mode = a.substr(std::strlen("--prime="));
            if (opt.prime_mode != "31" && opt.prime_mode != "61") {
                throw std::invalid_argument("unknown --prime=..., expected 31|61");
            }
        } else if (starts_with(a, "--min_pow=")) {
            opt.min_pow = static_cast<int>(parse_u64(a.substr(std::strlen("--min_pow="))));
        } else if (starts_with(a, "--max_pow=")) {
            opt.max_pow = static_cast<int>(parse_u64(a.substr(std::strlen("--max_pow="))));
        } else if (starts_with(a, "--repeats=")) {
            opt.repeats = static_cast<int>(parse_u64(a.substr(std::strlen("--repeats="))));
        } else if (starts_with(a, "--seed=")) {
            opt.seed = parse_u64(a.substr(std::strlen("--seed=")));
        } else if (starts_with(a, "--csv=")) {
            opt.csv = (parse_u64(a.substr(std::strlen("--csv="))) != 0);
        } else if (starts_with(a, "--verify=")) {
            opt.verify = (parse_u64(a.substr(std::strlen("--verify="))) != 0);
        } else {
            throw std::invalid_argument("unknown argument: " + a);
        }
    }

    if (opt.min_pow < 0 || opt.max_pow < opt.min_pow) {
        throw std::invalid_argument("invalid min_pow/max_pow");
    }
    if (opt.repeats <= 0) {
        throw std::invalid_argument("repeats must be >= 1");
    }
    return opt;
}

static void run_one_prime(u64 p, const std::string& prime_label, const Options& opt) {
    FpCtx F(p);

    if (opt.csv) {
        std::cout << "prime,n_pow,n,base_bits,baseline_avg_ms,kronecker_avg_ms,ratio,baseline_fp,kronecker_fp\n";
    } else {
        std::cout << "Prime " << prime_label << " (p=" << p << "), repeats=" << opt.repeats << "\n";
    }

    std::mt19937_64 rng(opt.seed ^ (p + 0x9e3779b97f4a7c15ULL));
    std::uniform_int_distribution<u64> dist_nonzero(1, p - 1);

    for (int k = opt.min_pow; k <= opt.max_pow; ++k) {
        const std::size_t n = (std::size_t)1ULL << (unsigned)k;

        std::vector<FpPoly> a, b;
        a.reserve(opt.repeats);
        b.reserve(opt.repeats);
        for (int rep = 0; rep < opt.repeats; ++rep) {
            a.push_back(random_poly(F, n, rng, dist_nonzero));
            b.push_back(random_poly(F, n, rng, dist_nonzero));
        }

        Timings base = bench_mul_baseline(a, b);
        unsigned base_bits = 0;
        Timings kron = bench_mul_kronecker(a, b, base_bits);

        if (opt.verify) {
            FpPoly check = a[0].mul(b[0]);
            FpPoly check_k = mul_kronecker(a[0], b[0]);
            if (!check.equals(check_k)) {
                std::cerr << "Verification failed at n=2^" << k << "\n";
                std::exit(1);
            }
        }

        if (opt.csv) {
            std::cout << prime_label << ","
                      << k << "," << n << ","
                      << base_bits << ","
                      << std::fixed << std::setprecision(3)
                      << base.avg_ms << ","
                      << kron.avg_ms << ","
                      << (kron.avg_ms / base.avg_ms) << ","
                      << std::hex << base.fingerprint << ","
                      << kron.fingerprint << std::dec
                      << "\n";
        } else {
            std::cout << "  n=2^" << k << " (" << n << "), base_bits=" << base_bits
                      << ": baseline avg=" << std::fixed << std::setprecision(3) << base.avg_ms << " ms, "
                      << "kron avg=" << kron.avg_ms << " ms "
                      << "(ratio=" << std::setprecision(3) << (kron.avg_ms / base.avg_ms) << ")"
                      << ", fp=0x" << std::hex << kron.fingerprint << std::dec
                      << "\n";
        }
    }
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);

        if (opt.prime_mode == "31") {
            run_one_prime(mersenne31(), "M31", opt);
        } else if (opt.prime_mode == "61") {
            run_one_prime(mersenne61(), "M61", opt);
        } else {
            throw std::logic_error("unreachable prime_mode");
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        print_help(argv[0]);
        return 1;
    }
}
