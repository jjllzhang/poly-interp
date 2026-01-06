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

// ---------------------------------- helpers ----------------------------------
static u64 mersenne31() { return (1ULL << 31) - 1; } // 2^31 - 1
static u64 mersenne61() { return (1ULL << 61) - 1; } // 2^61 - 1

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

// Lightweight polynomial fingerprint to keep compiler from optimizing results away.
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

static u64 vec_fingerprint(const std::vector<Fp>& v) {
    u64 h = 1469598103934665603ULL;
    h = mix64(h, static_cast<u64>(v.size()));
    for (std::size_t i = 0; i < v.size(); i += 1 + (v.size() / 7)) { // sample a few entries
        h = mix64(h, v[i].v);
    }
    return h;
}

// ---------------------------------- options ----------------------------------
struct Options {
    // prime selection: 31 / 61 / all / custom
    std::string prime_mode = "all";
    u64 custom_p = 0;

    int min_pow = 8;  // n = 2^k
    int max_pow = 14;

    int repeats = 3;  // per operation
    u64 seed = 1234567;

    bool csv = true;
};

static void print_help(const char* argv0) {
    std::cout
        << "Usage:\n"
        << "  " << argv0 << " [options]\n\n"
        << "Options:\n"
        << "  --prime=31|61|all        choose Mersenne prime (default: all)\n"
        << "  --p=<uint64>             custom prime modulus (no primality check)\n"
        << "  --min_pow=<k>            min k where n=2^k (default: 8)\n"
        << "  --max_pow=<k>            max k where n=2^k (default: 14)\n"
        << "  --repeats=<r>            repeats per operation (default: 3)\n"
        << "  --seed=<s>               RNG seed (default: 1234567)\n"
        << "  --csv=0|1                output CSV (default: 1)\n"
        << "  --help                   show this help\n\n"
        << "Notes:\n"
        << "  This benchmark isolates core primitives used by interpolation: subproduct tree build,\n"
        << "  polynomial mul/div/mod, derivative, multipoint eval, and batch inversion. Large n may\n"
        << "  still take noticeable time with current naive/NTT mix; adjust --max_pow accordingly.\n";
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

// ---------------------------------- benchmark ----------------------------------
struct Timings {
    double avg_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
    u64 fingerprint = 0;
};

// Generate a random polynomial of given length; forces leading coeff to be non-zero.
static FpPoly random_poly(const FpCtx& F, std::size_t len, std::mt19937_64& rng,
                          std::uniform_int_distribution<u64>& dist_nonzero) {
    std::vector<Fp> c(len);
    for (std::size_t i = 0; i < len; ++i) c[i] = Fp{dist_nonzero(rng)}; // in [1, p-1]
    c[0].v %= F.modulus(); // constant term still valid
    return FpPoly(F, std::move(c));
}

static Timings bench_tree_build(const FpCtx& F,
                                const std::vector<Fp>& xs,
                                int repeats) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        auto tree = FpPoly::SubproductTree::build(F, xs);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(tree.root()), (u64)rep + 1);
        sink ^= t.fingerprint;
    }
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_poly_mul(const std::vector<FpPoly>& a,
                              const std::vector<FpPoly>& b) {
    if (a.size() != b.size()) throw std::logic_error("bench_poly_mul: size mismatch");
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(a.size());

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

static Timings bench_divrem(const std::vector<FpPoly>& dividends,
                            const std::vector<FpPoly>& divisors) {
    if (dividends.size() != divisors.size()) throw std::logic_error("bench_divrem: size mismatch");
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(dividends.size());

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        auto qr = dividends[rep].divrem(divisors[rep]);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        u64 fp = mix64(poly_fingerprint(qr.first), poly_fingerprint(qr.second));
        t.fingerprint ^= mix64(fp, (u64)rep + 1);
        sink ^= t.fingerprint;
    }
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_mod(const std::vector<FpPoly>& dividends,
                         const std::vector<FpPoly>& divisors) {
    if (dividends.size() != divisors.size()) throw std::logic_error("bench_mod: size mismatch");
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(dividends.size());

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        FpPoly r = dividends[rep].mod(divisors[rep]);
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

static Timings bench_derivative(const std::vector<FpPoly>& polys) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(polys.size());

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        FpPoly d = polys[rep].derivative();
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(d), (u64)rep + 1);
        sink ^= t.fingerprint;
    }
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_multipoint_eval(const std::vector<FpPoly>& polys,
                                     const FpPoly::SubproductTree& tree) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(polys.size());

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        auto ys = polys[rep].multipoint_eval(tree);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(vec_fingerprint(ys), (u64)rep + 1);
        sink ^= t.fingerprint;
    }
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_batch_inv(const FpCtx& F,
                               std::vector<std::vector<Fp>> inputs) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const int repeats = static_cast<int>(inputs.size());

    for (int rep = 0; rep < repeats; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        F.batch_inv(inputs[rep]);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(vec_fingerprint(inputs[rep]), (u64)rep + 1);
        sink ^= t.fingerprint;
    }
    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static void print_row_csv(const std::string& prime_label,
                          const std::string& op,
                          int k,
                          std::size_t n,
                          const Timings& t) {
    const double ms_per_item = t.avg_ms / static_cast<double>(n);
    std::cout << prime_label << ","
              << op << ","
              << k << ","
              << n << ","
              << std::fixed << std::setprecision(3)
              << t.avg_ms << ","
              << t.min_ms << ","
              << t.max_ms << ","
              << std::setprecision(9) << ms_per_item << ","
              << std::hex << t.fingerprint << std::dec
              << "\n";
}

static void run_one_prime(u64 p, const std::string& prime_label, const Options& opt, bool print_header) {
    FpCtx F(p);

    if (opt.csv && print_header) {
        std::cout << "prime,op,n_pow,n,avg_ms,min_ms,max_ms,ms_per_item,fingerprint\n";
    }
    if (!opt.csv) {
        std::cout << "Prime " << prime_label << " (p=" << p << ")\n";
        std::cout << "n=2^k, k in [" << opt.min_pow << "," << opt.max_pow << "], repeats=" << opt.repeats << "\n";
    }

    std::mt19937_64 rng(opt.seed ^ (p + 0x9e3779b97f4a7c15ULL));
    std::uniform_int_distribution<u64> dist_nonzero(1, p - 1);

    for (int k = opt.min_pow; k <= opt.max_pow; ++k) {
        const std::size_t n = (std::size_t)1ULL << (unsigned)k;

        // Prepare points x_i = i for tree-based routines
        std::vector<Fp> xs(n);
        for (std::size_t i = 0; i < n; ++i) xs[i] = Fp{ static_cast<u64>(i) };

        // Inputs for polynomial ops
        std::vector<FpPoly> mul_a, mul_b;
        std::vector<FpPoly> div_a, div_b;
        std::vector<FpPoly> eval_polys;
        mul_a.reserve(opt.repeats);
        mul_b.reserve(opt.repeats);
        div_a.reserve(opt.repeats);
        div_b.reserve(opt.repeats);
        eval_polys.reserve(opt.repeats);

        for (int rep = 0; rep < opt.repeats; ++rep) {
            mul_a.push_back(random_poly(F, n, rng, dist_nonzero));
            mul_b.push_back(random_poly(F, n, rng, dist_nonzero));

            div_b.push_back(random_poly(F, n, rng, dist_nonzero));       // divisor degree n-1
            div_a.push_back(random_poly(F, 2 * n - 1, rng, dist_nonzero)); // dividend degree 2n-2

            eval_polys.push_back(random_poly(F, n, rng, dist_nonzero));  // degree n-1
        }

        // Batch inversion inputs (non-zero)
        std::vector<std::vector<Fp>> inv_inputs;
        inv_inputs.reserve(opt.repeats);
        for (int rep = 0; rep < opt.repeats; ++rep) {
            std::vector<Fp> v(n);
            for (std::size_t i = 0; i < n; ++i) v[i] = Fp{ dist_nonzero(rng) };
            inv_inputs.push_back(std::move(v));
        }

        // Build one tree for multipoint evaluation (not timed)
        auto eval_tree = FpPoly::SubproductTree::build(F, xs);

        Timings tree = bench_tree_build(F, xs, opt.repeats);
        Timings mul  = bench_poly_mul(mul_a, mul_b);
        Timings div  = bench_divrem(div_a, div_b);
        Timings mod  = bench_mod(div_a, div_b);
        Timings deriv = bench_derivative(eval_polys);
        Timings mpev = bench_multipoint_eval(eval_polys, eval_tree);
        Timings binv = bench_batch_inv(F, inv_inputs);

        if (opt.csv) {
            print_row_csv(prime_label, "subproduct_build", k, n, tree);
            print_row_csv(prime_label, "poly_mul", k, n, mul);
            print_row_csv(prime_label, "poly_divrem", k, n, div);
            print_row_csv(prime_label, "poly_mod", k, n, mod);
            print_row_csv(prime_label, "derivative", k, n, deriv);
            print_row_csv(prime_label, "multipoint_eval_tree", k, n, mpev);
            print_row_csv(prime_label, "batch_inv", k, n, binv);
        } else {
            auto pr = [&](const std::string& op, const Timings& t) {
                std::cout << "  [" << op << "] n=2^" << k << " (" << n << "): "
                          << std::fixed << std::setprecision(3)
                          << "avg=" << t.avg_ms << " ms, "
                          << "min=" << t.min_ms << ", max=" << t.max_ms
                          << ", ms/item=" << std::setprecision(9) << (t.avg_ms / (double)n)
                          << ", fp=0x" << std::hex << t.fingerprint << std::dec << "\n";
            };
            pr("subproduct_build", tree);
            pr("poly_mul", mul);
            pr("poly_divrem", div);
            pr("poly_mod", mod);
            pr("derivative", deriv);
            pr("multipoint_eval_tree", mpev);
            pr("batch_inv", binv);
        }
    }
}

int main(int argc, char** argv) {
    try {
        Options opt = parse_args(argc, argv);

        if (opt.prime_mode == "31") {
            run_one_prime(mersenne31(), "M31", opt, true);
        } else if (opt.prime_mode == "61") {
            run_one_prime(mersenne61(), "M61", opt, true);
        } else if (opt.prime_mode == "all") {
            run_one_prime(mersenne31(), "M31", opt, true);
            run_one_prime(mersenne61(), "M61", opt, false);
        } else if (opt.prime_mode == "custom") {
            run_one_prime(opt.custom_p, "custom", opt, true);
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
