#include <poly_interp/prime_field.hpp>
#include <poly_interp/fp_poly.hpp>


#include <chrono>
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
            tree = FpPoly::SubproductTree::build(F, xs);
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
                poly = FpPoly::interpolate_subproduct_tree(tree, ys);
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
