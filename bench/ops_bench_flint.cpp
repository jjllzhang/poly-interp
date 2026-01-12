#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>

#include <algorithm>
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

using u64 = std::uint64_t;

static u64 mersenne31() { return (1ULL << 31) - 1; }                 // 2^31 - 1
static u64 mersenne61() { return (1ULL << 61) - 1; }                 // 2^61 - 1

static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && std::memcmp(s.data(), prefix.data(), prefix.size()) == 0;
}

static u64 parse_u64(const std::string& s) {
    std::size_t pos = 0;
    unsigned long long v = std::stoull(s, &pos, 10);
    if (pos != s.size()) throw std::invalid_argument("invalid integer: " + s);
    return static_cast<u64>(v);
}

static inline u64 splitmix64(u64& x) {
    u64 z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static u64 mix64(u64 h, u64 x) {
    h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

// Lightweight fingerprint for nmod_poly to keep results live.
static u64 poly_fingerprint(const nmod_poly_t f) {
    slong deg = nmod_poly_degree(f);
    u64 h = 1469598103934665603ULL;
    h = mix64(h, static_cast<u64>(deg));
    if (deg < 0) return h;

    auto coeff = [&](slong i) -> u64 {
        if (i < 0) return 0;
        return nmod_poly_get_coeff_ui(f, static_cast<ulong>(i));
    };

    h = mix64(h, coeff(0));
    h = mix64(h, coeff(1));
    h = mix64(h, coeff(2));
    h = mix64(h, coeff(deg / 2));
    h = mix64(h, coeff(deg));
    return h;
}

static u64 vec_fingerprint(const std::vector<ulong>& v) {
    u64 h = 1469598103934665603ULL;
    h = mix64(h, static_cast<u64>(v.size()));
    if (v.empty()) return h;
    const std::size_t step = std::max<std::size_t>(1, v.size() / 7);
    for (std::size_t i = 0; i < v.size(); i += step) {
        h = mix64(h, v[i]);
    }
    return h;
}

struct Options {
    std::string prime_mode = "all"; // 31 | 61 | all | custom
    u64 custom_p = 0;
    int min_pow = 8;
    int max_pow = 14;
    int repeats = 3;
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
        << "  Benchmarks FLINT nmod_poly primitives: product_roots, mul, div/rem, mod,\n"
        << "  derivative, multipoint eval, and batch inversion on field elements.\n";
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

struct Timings {
    double avg_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
    u64 fingerprint = 0;
};

static void random_poly(nmod_poly_t f, ulong mod, slong len, u64& rng_state) {
    nmod_poly_zero(f);
    nmod_poly_fit_length(f, len);
    for (slong i = 0; i < len; ++i) {
        u64 c = splitmix64(rng_state) % mod;
        if (i == len - 1 && c == 0) c = 1; // force leading coeff non-zero
        nmod_poly_set_coeff_ui(f, static_cast<ulong>(i), static_cast<ulong>(c));
    }
    _nmod_poly_normalise(f);
}

static Timings bench_product_roots(const std::vector<ulong>& xs,
                                   ulong mod,
                                   int repeats) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;

    for (int rep = 0; rep < repeats; ++rep) {
        nmod_poly_t prod;
        nmod_poly_init2(prod, mod, xs.size() + 1);

        auto t0 = std::chrono::steady_clock::now();
        nmod_poly_product_roots_nmod_vec(prod, xs.data(), (slong)xs.size());
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(prod), (u64)rep + 1);
        sink ^= t.fingerprint;

        nmod_poly_clear(prod);
    }

    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_poly_mul(ulong mod, slong len, int repeats, u64& rng_state) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;

    for (int rep = 0; rep < repeats; ++rep) {
        nmod_poly_t a, b, r;
        nmod_poly_init2(a, mod, len);
        nmod_poly_init2(b, mod, len);
        nmod_poly_init(r, mod);
        random_poly(a, mod, len, rng_state);
        random_poly(b, mod, len, rng_state);

        auto t0 = std::chrono::steady_clock::now();
        nmod_poly_mul(r, a, b);
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(r), (u64)rep + 1);
        sink ^= t.fingerprint;

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r);
    }

    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_divrem(ulong mod, slong len_dividend, slong len_divisor, int repeats, u64& rng_state) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;

    for (int rep = 0; rep < repeats; ++rep) {
        nmod_poly_t a, b, q, r;
        nmod_poly_init2(a, mod, len_dividend);
        nmod_poly_init2(b, mod, len_divisor);
        nmod_poly_init(q, mod);
        nmod_poly_init(r, mod);
        random_poly(a, mod, len_dividend, rng_state);
        random_poly(b, mod, len_divisor, rng_state);

        auto t0 = std::chrono::steady_clock::now();
        nmod_poly_divrem(q, r, a, b);
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        u64 fp = mix64(poly_fingerprint(q), poly_fingerprint(r));
        t.fingerprint ^= mix64(fp, (u64)rep + 1);
        sink ^= t.fingerprint;

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(r);
    }

    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_mod(ulong mod, slong len_dividend, slong len_divisor, int repeats, u64& rng_state) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;

    for (int rep = 0; rep < repeats; ++rep) {
        nmod_poly_t a, b, r;
        nmod_poly_init2(a, mod, len_dividend);
        nmod_poly_init2(b, mod, len_divisor);
        nmod_poly_init(r, mod);
        random_poly(a, mod, len_dividend, rng_state);
        random_poly(b, mod, len_divisor, rng_state);

        auto t0 = std::chrono::steady_clock::now();
        nmod_poly_rem(r, a, b);
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(r), (u64)rep + 1);
        sink ^= t.fingerprint;

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(r);
    }

    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_derivative(ulong mod, slong len, int repeats, u64& rng_state) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;

    for (int rep = 0; rep < repeats; ++rep) {
        nmod_poly_t f, d;
        nmod_poly_init2(f, mod, len);
        nmod_poly_init(d, mod);
        random_poly(f, mod, len, rng_state);

        auto t0 = std::chrono::steady_clock::now();
        nmod_poly_derivative(d, f);
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(poly_fingerprint(d), (u64)rep + 1);
        sink ^= t.fingerprint;

        nmod_poly_clear(f);
        nmod_poly_clear(d);
    }

    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_multipoint_eval(const std::vector<ulong>& xs,
                                     ulong mod,
                                     slong len_poly,
                                     int repeats,
                                     u64& rng_state) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    const slong n = (slong)xs.size();

    for (int rep = 0; rep < repeats; ++rep) {
        nmod_poly_t f;
        nmod_poly_init2(f, mod, len_poly);
        random_poly(f, mod, len_poly, rng_state);

        std::vector<ulong> ys(xs.size(), 0);
        auto t0 = std::chrono::steady_clock::now();
        nmod_poly_evaluate_nmod_vec(ys.data(), f, xs.data(), n);
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(vec_fingerprint(ys), (u64)rep + 1);
        sink ^= t.fingerprint;

        nmod_poly_clear(f);
    }

    t.avg_ms /= static_cast<double>(repeats);
    if (sink == 0xdeadbeefULL) std::cerr << "sink\n";
    return t;
}

static Timings bench_batch_inv(const std::vector<ulong>& xs_template,
                               ulong mod,
                               int repeats,
                               u64& rng_state) {
    Timings t;
    t.min_ms = 1e300;
    volatile u64 sink = 0;
    nmod_t ctx;
    nmod_init(&ctx, mod);

    const std::size_t n = xs_template.size();
    std::vector<ulong> vals(n);
    std::vector<ulong> pref(n);

    for (int rep = 0; rep < repeats; ++rep) {
        // fresh random non-zero inputs
        for (std::size_t i = 0; i < n; ++i) {
            ulong v = (ulong)(splitmix64(rng_state) % mod);
            if (v == 0) v = 1;
            vals[i] = v;
        }

        auto t0 = std::chrono::steady_clock::now();
        pref[0] = vals[0];
        for (std::size_t i = 1; i < n; ++i) {
            pref[i] = nmod_mul(pref[i - 1], vals[i], ctx);
        }

        ulong inv_total = nmod_inv(pref[n - 1], ctx);
        for (std::size_t i = n; i-- > 0;) {
            ulong before = (i == 0) ? 1UL : pref[i - 1];
            ulong old = vals[i];
            vals[i] = nmod_mul(inv_total, before, ctx);
            inv_total = nmod_mul(inv_total, old, ctx);
        }
        auto t1 = std::chrono::steady_clock::now();

        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        t.avg_ms += ms;
        t.min_ms = std::min(t.min_ms, ms);
        t.max_ms = std::max(t.max_ms, ms);
        t.fingerprint ^= mix64(vec_fingerprint(vals), (u64)rep + 1);
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

static void run_one_prime(ulong p, const std::string& prime_label, const Options& opt, bool print_header) {
    if (opt.csv && print_header) {
        std::cout << "prime,op,n_pow,n,avg_ms,min_ms,max_ms,ms_per_item,fingerprint\n";
    }
    if (!opt.csv) {
        std::cout << "Prime " << prime_label << " (p=" << p << ")\n";
        std::cout << "n=2^k, k in [" << opt.min_pow << "," << opt.max_pow << "], repeats=" << opt.repeats << "\n";
    }

    // RNG
    u64 rng_state = opt.seed ^ (p + 0x9e3779b97f4a7c15ULL);

    for (int k = opt.min_pow; k <= opt.max_pow; ++k) {
        const std::size_t n = (std::size_t)1ULL << (unsigned)k;

        // xs = 0..n-1 (all < p for tested primes/sizes)
        std::vector<ulong> xs(n);
        for (std::size_t i = 0; i < n; ++i) xs[i] = (ulong)i % p;

        Timings prod = bench_product_roots(xs, p, opt.repeats);
        Timings mul = bench_poly_mul(p, (slong)n, opt.repeats, rng_state);
        Timings div = bench_divrem(p, (slong)(2 * n - 1), (slong)n, opt.repeats, rng_state);
        Timings rem = bench_mod(p, (slong)(2 * n - 1), (slong)n, opt.repeats, rng_state);
        Timings deriv = bench_derivative(p, (slong)n, opt.repeats, rng_state);
        Timings eval = bench_multipoint_eval(xs, p, (slong)n, opt.repeats, rng_state);
        Timings binv = bench_batch_inv(xs, p, opt.repeats, rng_state);

        if (opt.csv) {
            // Align op naming with project benchmarks for plotting/compare.
            print_row_csv(prime_label, "subproduct_build", k, n, prod);
            print_row_csv(prime_label, "poly_mul", k, n, mul);
            print_row_csv(prime_label, "poly_divrem", k, n, div);
            print_row_csv(prime_label, "poly_mod", k, n, rem);
            print_row_csv(prime_label, "derivative", k, n, deriv);
            print_row_csv(prime_label, "multipoint_eval", k, n, eval);
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
            pr("subproduct_build", prod);
            pr("poly_mul", mul);
            pr("poly_divrem", div);
            pr("poly_mod", rem);
            pr("derivative", deriv);
            pr("multipoint_eval", eval);
            pr("batch_inv", binv);
        }
    }
}

int main(int argc, char** argv) {
    Options opt = parse_args(argc, argv);

    if (opt.prime_mode == "31") {
        run_one_prime((ulong)mersenne31(), "M31", opt, true);
    } else if (opt.prime_mode == "61") {
        run_one_prime((ulong)mersenne61(), "M61", opt, true);
    } else if (opt.prime_mode == "all") {
        run_one_prime((ulong)mersenne31(), "M31", opt, true);
        run_one_prime((ulong)mersenne61(), "M61", opt, false);
    } else if (opt.prime_mode == "custom") {
        run_one_prime((ulong)opt.custom_p, "custom", opt, true);
    } else {
        std::cerr << "unreachable prime_mode\n";
        return 1;
    }

    return 0;
}
