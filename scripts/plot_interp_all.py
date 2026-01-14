"""
Plot interpolation benchmarks from multiple sources (project, naive, FLINT, NTL) on one set of charts.

Usage:
  python3 scripts/plot_interp_all.py \
      --bench data/bench.csv \
      --bench-naive data/bench_naive.csv \
      --bench-flint data/bench_flint.csv \
      --bench-ntl data/bench_NTL.csv \
      --output plots/interp_all
"""

import argparse
import os
import re
from typing import Optional, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def _load_csv(path: str) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path)
    except FileNotFoundError:
        print(f"Warning: file not found '{path}', skipping.")
        return None


def _normalize_interp_df(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """
    Normalize various interp CSV layouts into a common schema:
    columns: source, prime, n, build_time_ms, interp_time_ms, total_time_ms.
    """
    if df is None:
        return pd.DataFrame()

    df = df.copy()
    if "prime" in df.columns:
        df = df[df["prime"] != "prime"].copy()
        df["prime"] = df["prime"].astype(str).str.upper()
    else:
        df["prime"] = "UNKNOWN"

    numeric_cols = [
        "n_pow",
        "n",
        "degree",
        "build_tree_ms",
        "build_ms",
        "interp_pre_ms",
        "interp_core_avg_ms",
        "interp_core_min_ms",
        "interp_core_max_ms",
        "interp_avg_ms",
        "interp_min_ms",
        "interp_max_ms",
        "interp_ms_per_point",
        "interp_ms",
        "differentiate_ms",
        "evaluate_diff_ms",
        "interpolation_ms",
        "total_ms",
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if "n" not in df.columns:
        if "degree" in df.columns:
            df["n"] = df["degree"] + 1
        elif "n_pow" in df.columns:
            df["n"] = 2 ** df["n_pow"]

    build_col = "build_tree_ms" if "build_tree_ms" in df.columns else ("build_ms" if "build_ms" in df.columns else None)
    interp_col = "interp_avg_ms" if "interp_avg_ms" in df.columns else ("interp_ms" if "interp_ms" in df.columns else None)
    if build_col:
        df["build_time_ms"] = df[build_col]
    if interp_col:
        df["interp_time_ms"] = df[interp_col]
    if "total_ms" in df.columns:
        df["total_time_ms"] = df["total_ms"]
    elif build_col and interp_col:
        df["total_time_ms"] = df["build_time_ms"] + df["interp_time_ms"]

    df["source"] = source
    required = ["n", "build_time_ms", "interp_time_ms"]
    df = df.dropna(subset=required).copy()
    if "total_time_ms" not in df.columns:
        df["total_time_ms"] = df["build_time_ms"] + df["interp_time_ms"]

    keep_cols = ["source", "prime", "n", "build_time_ms", "interp_time_ms", "total_time_ms"]
    extra_cols = [c for c in ["interp_ms_per_point"] if c in df.columns]
    missing = [c for c in keep_cols if c not in df.columns]
    if missing:
        print(f"Warning: skipping source '{source}' due to missing columns: {', '.join(missing)}")
        return pd.DataFrame()
    return df[keep_cols + extra_cols]


def plot_metric(df: pd.DataFrame, metric: str, ylabel: str, output_dir: str, logy: bool) -> None:
    if metric not in df.columns:
        return
    for prime, g in df.groupby("prime"):
        plt.figure(figsize=(10, 6))
        for src, gg in g.groupby("source"):
            gg = gg.sort_values("n")
            plt.plot(gg["n"], gg[metric], marker="o", label=src)
        plt.xscale("log", base=2)
        if logy:
            plt.yscale("log")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.xlabel("n (number of points)")
        plt.ylabel(ylabel)
        plt.title(f"{metric} vs n (prime={prime})")
        plt.legend()
        plt.tight_layout()

        safe_prime = re.sub(r"[^A-Za-z0-9_-]+", "_", prime)
        fname = f"{metric}_{safe_prime}.png"
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, fname)
        plt.savefig(path)
        print(f"Saved: {path}")
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot interpolation benchmarks from project/naive/FLINT/NTL CSVs.")
    parser.add_argument("--bench", default="data/bench.csv", help="Project bench CSV (subproduct-tree).")
    parser.add_argument("--bench-naive", default="data/bench_naive.csv", help="Project naive Lagrange bench CSV.")
    parser.add_argument("--bench-flint", default="data/bench_flint.csv", help="FLINT bench CSV.")
    parser.add_argument("--bench-ntl", default="data/bench_NTL.csv", help="NTL bench CSV.")
    parser.add_argument("--output-dir", default="plots/interp_all", help="Output directory for figures.")
    parser.add_argument("--no-logy", action="store_true", help="Disable log scale on Y axis.")
    args = parser.parse_args()

    sources: List[Tuple[str, str]] = [
        ("Subproduct-tree(Ours)", args.bench),
        ("Lagrange", args.bench_naive),
        ("flint", args.bench_flint),
        ("ntl", args.bench_ntl),
    ]

    frames = []
    for label, path in sources:
        df_raw = _load_csv(path)
        if df_raw is None:
            continue
        df_norm = _normalize_interp_df(df_raw, label)
        if df_norm.empty:
            print(f"Warning: no valid rows for '{label}' ({path}), skipping.")
            continue
        frames.append(df_norm)

    if not frames:
        print("Error: no input data to plot.")
        return

    df_all = pd.concat(frames, ignore_index=True)
    df_all = df_all.sort_values(["prime", "n", "source"])

    logy = not args.no_logy
    metrics = [
        ("build_time_ms", "build time (ms)"),
        ("interp_time_ms", "interp time (ms)"),
        ("total_time_ms", "total time (ms)"),
    ]
    if "interp_ms_per_point" in df_all.columns:
        metrics.append(("interp_ms_per_point", "interp time per point (ms)"))

    for metric, ylabel in metrics:
        plot_metric(df_all, metric, ylabel, args.output_dir, logy)


if __name__ == "__main__":
    main()
