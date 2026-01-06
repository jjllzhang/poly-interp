#include <flint/flint.h>
#include <flint/nmod_poly.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#if __has_include(<filesystem>)
  #include <filesystem>
  namespace fs = std::filesystem;
#endif

static volatile ulong g_sink = 0;

static inline uint64_t splitmix64(uint64_t &x) {
  uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static void fill_xs(std::vector<ulong> &xs, ulong mod) {
  for (size_t i = 0; i < xs.size(); ++i) xs[i] = (ulong)((i + 1) % mod);
}

static void fill_ys(std::vector<ulong> &ys, ulong mod, uint64_t &seed) {
  for (size_t i = 0; i < ys.size(); ++i) ys[i] = (ulong)(splitmix64(seed) % mod);
}

static bool file_empty_or_missing(const std::string& path) {
#if __has_include(<filesystem>)
  std::error_code ec;
  if (!fs::exists(path, ec)) return true;
  auto sz = fs::file_size(path, ec);
  return ec ? true : (sz == 0);
#else
  std::ifstream in(path, std::ios::binary);
  return !in.good() || (in.peek() == std::ifstream::traits_type::eof());
#endif
}

int main(int argc, char **argv) {
#if FLINT_BITS < 64
  std::cerr << "FLINT_BITS < 64: cannot benchmark M61 on this build.\n";
#endif

  const ulong M31 = 2147483647UL;                 // 2^31 - 1
  const ulong M61 = 2305843009213693951UL;        // 2^61 - 1

  std::string mod_arg = "all";
  std::string csv_path;

  // Args: [m31|m61|all] [--csv out.csv]
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "m31" || a == "m61" || a == "all") mod_arg = a;
    else if (a == "--csv" && i + 1 < argc) csv_path = argv[++i];
  }

  std::vector<std::pair<std::string, ulong>> mods;
  if (mod_arg == "m31") mods.push_back({"m31", M31});
  else if (mod_arg == "m61") mods.push_back({"m61", M61});
  else if (mod_arg == "all") {
    mods.push_back({"m31", M31});
    mods.push_back({"m61", M61});
  } else {
    std::cerr << "Unknown modulus arg: " << mod_arg << "\n";
    return 1;
  }

  int k_min = 10, k_max = 20;
  int trials = 3;
  int warmup = 1;

  std::ofstream csv;
  if (!csv_path.empty()) {
    bool need_header = file_empty_or_missing(csv_path);
    csv.open(csv_path, std::ios::app);
    if (!csv) {
      std::cerr << "Failed to open CSV: " << csv_path << "\n";
      return 1;
    }
    if (need_header) {
      csv << "prime,label,n_pow,n,build_tree_ms,interp_avg_ms,interp_min_ms,interp_max_ms,interp_ms_per_point,fingerprint\n";
    }
  }

  for (auto [mod_name, mod] : mods) {
    const std::string prime_label = mod_name;
    std::cout << "Prime " << mod_name << " (p=" << mod << ")\n";
    std::cout << "n=2^k, k in [" << k_min << "," << k_max << "], trials=" << trials
              << " (avg over distinct ys)\n";

    uint64_t seed = 123456789ULL ^ (uint64_t)mod;

    for (int k = k_min; k <= k_max; ++k) {
      const slong n = (slong)1 << k;

      std::vector<ulong> xs((size_t)n), ys((size_t)n);
      fill_xs(xs, mod);

      nmod_poly_t poly;
      nmod_poly_init(poly, mod);

      // warmup
      for (int w = 0; w < warmup; ++w) {
        fill_ys(ys, mod, seed);
        nmod_poly_interpolate_nmod_vec_fast(poly, xs.data(), ys.data(), n);
        g_sink ^= nmod_poly_get_coeff_ui(poly, 0);
      }

      double sum_interp_ms = 0.0;
      double min_interp_ms = 1e300;
      double max_interp_ms = 0.0;
      double sum_tree_ms = 0.0;
      ulong fp_acc = 0;

      for (int t = 0; t < trials; ++t) {
        fill_ys(ys, mod, seed);

        // subproduct tree build timing via internal API (includes alloc/free)
        {
          auto t0 = std::chrono::steady_clock::now();
          nmod_t ctx;
          nmod_init(&ctx, mod);
          mp_ptr* tree = _nmod_poly_tree_alloc(n);
          _nmod_poly_tree_build(tree, xs.data(), n, ctx);
          _nmod_poly_tree_free(tree, n);
          auto t1 = std::chrono::steady_clock::now();
          sum_tree_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
        }

        auto st = std::chrono::steady_clock::now();
        nmod_poly_interpolate_nmod_vec_fast(poly, xs.data(), ys.data(), n);
        auto ed = std::chrono::steady_clock::now();

        // digest to avoid DCE
        ulong d = 0;
        slong len = nmod_poly_length(poly);
        for (slong i = 0; i < std::min<slong>(len, 8); ++i) {
          d ^= nmod_poly_get_coeff_ui(poly, (ulong)i) + 0x9e3779b9U;
        }
        g_sink ^= d;
        fp_acc ^= (d + 0x9e3779b97f4a7c15ULL + (fp_acc << 6) + (fp_acc >> 2));

        double ms = std::chrono::duration<double, std::milli>(ed - st).count();
        sum_interp_ms += ms;
        min_interp_ms = std::min(min_interp_ms, ms);
        max_interp_ms = std::max(max_interp_ms, ms);
      }

      double avg_interp_ms = sum_interp_ms / (double)trials;
      double avg_tree_ms = sum_tree_ms / (double)trials;
      double ms_per_point = avg_interp_ms / (double)n;

      if (csv) {
        csv
            << prime_label << "," << prime_label << ","
            << k << "," << n << ","
            << std::fixed << std::setprecision(3)
            << avg_tree_ms << ","
            << avg_interp_ms << ","
            << min_interp_ms << ","
            << max_interp_ms << ","
            << std::setprecision(9) << ms_per_point << ","
            << std::hex << fp_acc << std::dec
            << "\n";
        csv.flush();
      }

      if (csv_path.empty()) {
        std::cout
            << "n=2^" << k << " (" << n << "): "
            << "build=" << std::fixed << std::setprecision(3) << avg_tree_ms << " ms, "
            << "interp(avg)=" << avg_interp_ms << " ms "
            << "(min=" << min_interp_ms << ", max=" << max_interp_ms << "), "
            << "ms/pt=" << std::setprecision(9) << ms_per_point
            << ", fp=0x" << std::hex << fp_acc << std::dec
            << "\n";
      }

      nmod_poly_clear(poly);
    }
  }

  std::cerr << "sink=" << g_sink << "\n";
  return 0;
}
