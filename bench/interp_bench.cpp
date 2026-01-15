#include <poly_interp/prime_field.hpp>
#include <poly_interp/fp_poly.hpp>

#ifdef PI_HAVE_GMP
#include <gmp.h>
#endif

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

enum class InterpMethod {
    SubproductTree,
    NaiveLagrange
};

static u64 mersenne31() { return (1ULL << 31) - 1; }                 // 2^31 - 1
static u64 mersenne61() { return (1ULL << 61) - 1; }                 // 2^61 - 1

static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && std::memcmp(s.data(), prefix.data(), prefix.size()) == 0;
}

static u64 parse_u64(const std::string& s) {
    // 允许十进制
    std::size_t pos = 0;
    unsigned long long v = std::stoull(s, &pos, 10);
    if (pos != s.size()) throw std::invalid_argument("invalid integer: " + s);
    return static_cast<u64>(v);
}

struct Options {
    // prime selection:
    // prime_mode = "31" / "61" / "all" / "custom"
    std::string prime_mode = "61";
    u64 custom_p = 0;
    InterpMethod method = InterpMethod::SubproductTree;

    int min_pow = 10;
    int max_pow = 20;

    int repeats = 3;        // 插值重复次数（同一棵树、多组 ys），取平均
    int verify_samples = 0; // 抽样验证点数（0=不验证）

    u64 seed = 1234567;

    bool csv = true;        // 输出 CSV
};

static void print_help(const char* argv0) {
    std::cout
        << "Usage:\n"
        << "  " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --prime=31|61|all        choose Mersenne prime (default: 61)\n"
        << "  --p=<uint64>             custom prime modulus (no primality check)\n"
        << "  --method=tree|naive      interpolation method (default: tree)\n"
        << "  --min_pow=<k>            min k where n=2^k (default: 10)\n"
        << "  --max_pow=<k>            max k where n=2^k (default: 20)\n"
        << "  --repeats=<r>            repeats for interpolation per n (default: 3)\n"
        << "  --verify=<t>             verify t sampled points (default: 0)\n"
        << "  --seed=<s>               RNG seed (default: 1234567)\n"
        << "  --csv=0|1                output CSV (default: 1)\n"
        << "  --help                   show this help\n\n"
        << "Notes:\n"
        << "  --method=naive runs O(n^2) Lagrange interpolation; use smaller n (e.g. k<=14) to avoid long runtimes.\n";
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
            if (opt.prime_mode != "31" && opt.prime_mode != "61" && opt.prime_mode != "all") {
                throw std::invalid_argument("unknown --prime=..., expected 31|61|all");
            }
        } else if (starts_with(a, "--p=")) {
            opt.prime_mode = "custom";
            opt.custom_p = parse_u64(a.substr(std::strlen("--p=")));
        } else if (starts_with(a, "--method=")) {
            std::string m = a.substr(std::strlen("--method="));
            if (m == "tree") {
                opt.method = InterpMethod::SubproductTree;
            } else if (m == "naive") {
                opt.method = InterpMethod::NaiveLagrange;
            } else {
                throw std::invalid_argument("unknown --method=..., expected tree|naive");
            }
        } else if (starts_with(a, "--min_pow=")) {
            opt.min_pow = static_cast<int>(parse_u64(a.substr(std::strlen("--min_pow="))));
        } else if (starts_with(a, "--max_pow=")) {
            opt.max_pow = static_cast<int>(parse_u64(a.substr(std::strlen("--max_pow="))));
        } else if (starts_with(a, "--repeats=")) {
            opt.repeats = static_cast<int>(parse_u64(a.substr(std::strlen("--repeats="))));
        } else if (starts_with(a, "--verify=")) {
            opt.verify_samples = static_cast<int>(parse_u64(a.substr(std::strlen("--verify="))));
        } else if (starts_with(a, "--seed=")) {
            opt.seed = parse_u64(a.substr(std::strlen("--seed=")));
        } else if (starts_with(a, "--csv=")) {
            opt.csv = (parse_u64(a.substr(std::strlen("--csv="))) != 0);
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
    if (opt.verify_samples < 0) {
        throw std::invalid_argument("verify must be >= 0");
    }
    return opt;
}

// 生成一个轻量 fingerprint，避免编译器把结果优化掉，同时不引入 O(n) 额外开销
static u64 poly_fingerprint(const FpPoly& f) {
    // FNV-1a
    u64 h = 1469598103934665603ULL;
    auto mix = [&](u64 x) {
        h ^= x;
        h *= 1099511628211ULL;
    };

    int deg = f.degree();
    mix(static_cast<u64>(deg));

    if (deg < 0) return h;

    // 取常数、一次、二次、middle、leading（存在则取）
    auto pick = [&](std::size_t i) {
        Fp ci = f.coeff(i);
        mix(ci.v);
    };

    pick(0);
    if (deg >= 1) pick(1);
    if (deg >= 2) pick(2);

    std::size_t mid = static_cast<std::size_t>(deg) / 2;
    pick(mid);

    // leading
    mix(f.leading_coeff().v);

    return h;
}

#ifdef PI_HAVE_GMP
using u128 = unsigned __int128;

static unsigned choose_base_bits(u64 modulus, std::size_t min_len) {
    // max coeff <= min_len * (p-1)^2. Add a couple safety bits.
    double lp = std::log2(static_cast<double>(modulus));
    double bound = 2.0 * lp + std::log2(static_cast<double>(std::max<std::size_t>(min_len, 1))) + 2.0;
    unsigned bits = static_cast<unsigned>(std::ceil(bound));
    if (bits < 16) bits = 16;
    return bits;
}

static inline mp_limb_t mask_bits(unsigned bits) {
    if (bits == 0) return 0;
    if (bits >= GMP_NUMB_BITS) return ~(mp_limb_t)0;
    return ((mp_limb_t)1 << bits) - 1;
}

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
    const unsigned limbs_per_coeff = (base_bits % limb_bits == 0) ? (base_bits / limb_bits) : 0;
    if (limbs_per_coeff > 0) {
        const std::size_t limb_count = coeffs.size() * (std::size_t)limbs_per_coeff;
        out.assign(limb_count, 0);
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            out[i * (std::size_t)limbs_per_coeff] = (mp_limb_t)coeffs[i].v;
        }
        mp_size_t used = (mp_size_t)limb_count;
        while (used > 0 && out[(std::size_t)used - 1] == 0) --used;
        if (used == 0) used = 1;
        return used;
    }

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

    const unsigned limbs_per_coeff = (base_bits % limb_bits == 0) ? (base_bits / limb_bits) : 0;
    if (limbs_per_coeff > 0) {
        for (std::size_t i = 0; i < out.size(); ++i) {
            std::size_t idx = i * (std::size_t)limbs_per_coeff;
            if (idx >= limb_count_s) {
                out[i] = F.zero();
            } else {
                out[i] = F.from_uint((u64)limbs[idx] % F.modulus());
            }
        }
        return;
    }

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

static FpPoly mul_kronecker(const FpPoly& a, const FpPoly& b) {
    const FpCtx& F = a.ctx();
    if (a.modulus() != b.modulus()) {
        throw std::invalid_argument("mul_kronecker: modulus mismatch");
    }

    if (a.is_zero() || b.is_zero()) return FpPoly(F);

    const std::size_t n = a.coeffs().size();
    const std::size_t m = b.coeffs().size();
    const std::size_t min_len = std::min(n, m);
    const unsigned base_bits = choose_base_bits(F.modulus(), min_len);

    static thread_local MpnScratch S;
    static thread_local PackScratch PS;

    const mp_size_t limbs_a = pack_poly_bits(S.A, a.coeffs(), base_bits);
    const mp_size_t limbs_b = pack_poly_bits(S.B, b.coeffs(), base_bits);

    const mp_size_t limbs_res = limbs_a + limbs_b;
    S.C.resize((std::size_t)limbs_res);

    const bool same_operands = (&a == &b) || (a.coeffs().data() == b.coeffs().data() && a.coeffs().size() == b.coeffs().size());
    if (same_operands) {
        mpn_sqr(S.C.data(), S.A.data(), limbs_a);
    } else if (limbs_a == limbs_b) {
        mpn_mul_n(S.C.data(), S.A.data(), S.B.data(), limbs_a);
    } else if (limbs_a > limbs_b) {
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
#else
static FpPoly mul_kronecker(const FpPoly& a, const FpPoly& b) {
    // Fallback when GMP is unavailable.
    return a.mul(b);
}
#endif

static FpPoly::SubproductTree build_tree_kronecker(const FpCtx& ctx, const std::vector<Fp>& xs) {
    FpPoly::SubproductTree T(ctx);
    T.points = xs;
    for (auto& x : T.points) x.v %= ctx.modulus();

    if (T.points.empty()) return T;

    std::vector<FpPoly> level0;
    level0.reserve(T.points.size());
    for (const auto& xi : T.points) {
        FpPoly leaf(ctx);
        leaf.coeffs_mut().reserve(2);
        leaf.coeffs_mut().push_back(ctx.neg(xi));
        leaf.coeffs_mut().push_back(ctx.one());
        leaf.trim();
        level0.push_back(std::move(leaf));
    }
    T.levels.push_back(std::move(level0));

    while (T.levels.back().size() > 1) {
        const auto& prev = T.levels.back();
        std::vector<FpPoly> nxt;
        nxt.reserve((prev.size() + 1) / 2);

        for (std::size_t i = 0; i < prev.size(); i += 2) {
            if (i + 1 < prev.size()) {
                nxt.push_back(mul_kronecker(prev[i], prev[i + 1]));
            } else {
                nxt.push_back(prev[i]); // carry
            }
        }
        T.levels.push_back(std::move(nxt));
    }

    return T;
}

static FpPoly interpolate_subproduct_tree_kronecker(const FpPoly::SubproductTree& tree,
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

    const auto& inv_dvals = tree.inv_derivative_vals();

    std::vector<Fp> a(n, F.zero());
    for (std::size_t i = 0; i < n; ++i) {
        Fp yi = ys[i];
        yi.v %= F.modulus();
        a[i] = F.mul(yi, inv_dvals[i]);
    }

    std::vector<std::vector<FpPoly>> Flevels;
    Flevels.reserve(tree.n_levels());

    {
        std::vector<FpPoly> L0;
        L0.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            FpPoly leaf(F);
            if (a[i].v != 0) {
                leaf.coeffs_mut().push_back(a[i]);
            }
            leaf.trim();
            L0.push_back(std::move(leaf));
        }
        Flevels.push_back(std::move(L0));
    }

    for (std::size_t level = 1; level < tree.levels.size(); ++level) {
        const auto& prevF = Flevels[level - 1];
        const auto& prevM = tree.levels[level - 1];

        if (prevF.size() != prevM.size()) {
            throw std::invalid_argument("interpolate_subproduct_tree: tree malformed (prevF size != prevM size)");
        }

        std::vector<FpPoly> cur;
        cur.reserve((prevF.size() + 1) / 2);

        for (std::size_t i = 0; i < prevF.size(); i += 2) {
            if (i + 1 < prevF.size()) {
                FpPoly t1 = mul_kronecker(prevF[i], prevM[i + 1]);
                FpPoly t2 = mul_kronecker(prevF[i + 1], prevM[i]);
                FpPoly sum = t1.add(t2);
                sum.trim();
                cur.push_back(std::move(sum));
            } else {
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

static void run_one_prime(u64 p, const std::string& prime_label, const Options& opt) {
    FpCtx F(p);
    const bool method_naive = (opt.method == InterpMethod::NaiveLagrange);
    const std::string label = method_naive ? (prime_label + "-naive") : prime_label;

    // 输出表头
    if (opt.csv) {
        std::cout
            << "prime,label,n_pow,n,"
            << "build_tree_ms,"
            << "interp_pre_ms,"                 // precompute inv(dG(x_i)) once
            << "interp_core_avg_ms,"            // pure interpolation (reuse precompute)
            << "interp_core_min_ms,"
            << "interp_core_max_ms,"
            << "interp_avg_ms,"                 // total = pre + core
            << "interp_min_ms,"
            << "interp_max_ms,"
            << "interp_ms_per_point,"
            << "fingerprint\n";
    } else {
        std::cout << "Prime " << prime_label << " (p=" << p << "), method="
                  << (method_naive ? "naive-lagrange" : "subproduct-tree") << "\n";
        std::cout << "n=2^k, k in [" << opt.min_pow << "," << opt.max_pow
                  << "], repeats=" << opt.repeats << "\n";
    }

    std::mt19937_64 rng(opt.seed ^ (p + 0x9e3779b97f4a7c15ULL));
    std::uniform_int_distribution<u64> dist(0, p - 1);

    // 用于防止优化
    volatile u64 global_sink = 0;

    for (int k = opt.min_pow; k <= opt.max_pow; ++k) {
        const std::size_t n = (std::size_t)1ULL << (unsigned)k;

        // 生成点：x_i = i（保证互异，且 i < p 对 M31/M61 且 n<=2^20 成立）
        std::vector<Fp> xs(n);
        for (std::size_t i = 0; i < n; ++i) {
            xs[i] = Fp{ static_cast<u64>(i) };
        }

        // 1) build tree timing（naive 跳过）
        double build_ms = 0.0;
        FpPoly::SubproductTree tree;
        if (!method_naive) {
            auto t0 = std::chrono::steady_clock::now();
            tree = build_tree_kronecker(F, xs);
            auto t1 = std::chrono::steady_clock::now();
            build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        }

        // 2) interpolation timing
        double precompute_ms = 0.0;  // measured once per n
        bool precomputed = false;
        double core_sum_ms = 0.0;
        double core_min_ms = 1e300;
        double core_max_ms = 0.0;
        u64 fp_acc = 0;

        for (int rep = 0; rep < opt.repeats; ++rep) {
            // 生成 ys（不计入插值计时）
            std::vector<Fp> ys(n);
            for (std::size_t i = 0; i < n; ++i) {
                ys[i] = Fp{ dist(rng) }; // already in [0,p)
            }

            FpPoly poly(F);

            double core_ms = 0.0;
            if (!method_naive) {
                double pre_ms_this = 0.0;
                if (!precomputed) {
                    auto tpre0 = std::chrono::steady_clock::now();
                    tree.inv_derivative_vals();
                    auto tpre1 = std::chrono::steady_clock::now();
                    pre_ms_this = std::chrono::duration<double, std::milli>(tpre1 - tpre0).count();
                    precompute_ms = pre_ms_this;
                    precomputed = true;
                }
                auto tcore0 = std::chrono::steady_clock::now();
                poly = interpolate_subproduct_tree_kronecker(tree, ys);
                auto tcore1 = std::chrono::steady_clock::now();
                core_ms = std::chrono::duration<double, std::milli>(tcore1 - tcore0).count();
            } else {
                auto tcore0 = std::chrono::steady_clock::now();
                poly = FpPoly::interpolate_lagrange_naive(F, xs, ys);
                auto tcore1 = std::chrono::steady_clock::now();
                core_ms = std::chrono::duration<double, std::milli>(tcore1 - tcore0).count();
            }

            core_sum_ms += core_ms;
            core_min_ms = std::min(core_min_ms, core_ms);
            core_max_ms = std::max(core_max_ms, core_ms);

            // 抽样验证（不计入 interp_ms）
            if (opt.verify_samples > 0) {
                for (int t = 0; t < opt.verify_samples; ++t) {
                    // deterministic sample index
                    std::size_t idx = (std::size_t)((u64)t * 1315423911ULL + (u64)rep * 2654435761ULL) & (n - 1);
                    Fp v = poly.eval(xs[idx]);
                    if (v != ys[idx]) {
                        std::cerr << "Verification failed at n=2^" << k
                                  << ", rep=" << rep << ", idx=" << idx << "\n";
                        std::cerr << "  got=" << v.v << ", expected=" << ys[idx].v << "\n";
                        std::exit(1);
                    }
                }
            }

            // fingerprint（O(1)）
            u64 h = poly_fingerprint(poly);
            fp_acc ^= (h + 0x9e3779b97f4a7c15ULL + (fp_acc << 6) + (fp_acc >> 2));
            global_sink ^= fp_acc;
        }

        const double core_avg_ms = core_sum_ms / (double)opt.repeats;
        const double total_avg_ms = precompute_ms + core_avg_ms; // 预计算一次 + 一次核心插值的平均时间
        const double total_min_ms = precompute_ms + core_min_ms;
        const double total_max_ms = precompute_ms + core_max_ms;
        const double ms_per_point = total_avg_ms / (double)n;

        if (opt.csv) {
            std::cout
                << prime_label << "," << label << ","
                << k << "," << n << ","
                << std::fixed << std::setprecision(3)
                << build_ms << ","
                << precompute_ms << ","
                << core_avg_ms << ","
                << core_min_ms << ","
                << core_max_ms << ","
                << total_avg_ms << ","
                << total_min_ms << ","
                << total_max_ms << ","
                << std::setprecision(9) << ms_per_point << ","
                << std::hex << fp_acc << std::dec
                << "\n";
        } else {
            std::cout
                << "n=2^" << k << " (" << n << "): "
                << "build=" << std::fixed << std::setprecision(3) << build_ms << " ms, "
                << "interp_pre=" << precompute_ms << " ms, "
                << "interp_core(avg)=" << core_avg_ms << " ms "
                << "(min=" << core_min_ms << ", max=" << core_max_ms << "), "
                << "interp_total(avg)=" << total_avg_ms << " ms "
                << "(min=" << total_min_ms << ", max=" << total_max_ms << "), "
                << "ms/pt=" << std::setprecision(9) << ms_per_point
                << ", fp=0x" << std::hex << fp_acc << std::dec
                << "\n";
        }

        // 防止编译器“过分聪明”
        if (global_sink == 0xdeadbeefULL) std::cerr << "sink\n";
    }
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);

        if (opt.prime_mode == "31") {
            run_one_prime(mersenne31(), "M31", opt);
        } else if (opt.prime_mode == "61") {
            run_one_prime(mersenne61(), "M61", opt);
        } else if (opt.prime_mode == "all") {
            run_one_prime(mersenne31(), "M31", opt);
            run_one_prime(mersenne61(), "M61", opt);
        } else if (opt.prime_mode == "custom") {
            run_one_prime(opt.custom_p, "custom", opt);
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
