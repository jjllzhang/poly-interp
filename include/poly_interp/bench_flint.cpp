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
      csv << "lib,mod,k,n,avg_us\n";
    }
  }

  for (auto [mod_name, mod] : mods) {
    std::cout << "FLINT nmod_poly_interpolate_nmod_vec_fast  mod=" << mod
              << " (" << mod_name << ")  k=[" << k_min << "," << k_max << "]"
              << "  trials=" << trials << " (avg over distinct ys)\n";

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

      double sum_us = 0.0;

      for (int t = 0; t < trials; ++t) {
        fill_ys(ys, mod, seed);

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

        double us = std::chrono::duration<double, std::micro>(ed - st).count();
        sum_us += us;
      }

      double avg = sum_us / (double)trials;
      std::cout << "n=" << n << ", avg_us=" << avg << "\n";
      if (csv) {
        csv << "flint," << mod << "," << k << "," << n << "," << avg << "\n";
        csv.flush();
      }

      nmod_poly_clear(poly);
    }
  }

  std::cerr << "sink=" << g_sink << "\n";
  return 0;
}
