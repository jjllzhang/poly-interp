import argparse
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")  # no-GUI backend
import matplotlib.pyplot as plt
from typing import Optional


def _load_csv(csv_path: str) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: File '{csv_path}' not found.")
        return None


def plot_interp_csv(csv_path: str, logy: bool = True, output_dir: str = "plots") -> None:
    """
    Visualize interp_bench CSV and save figures to disk.
    """
    os.makedirs(output_dir, exist_ok=True)

    df = _load_csv(csv_path)
    if df is None:
        return

    # Drop repeated header rows if any
    if "prime" in df.columns:
        df = df[df["prime"] != "prime"].copy()

    num_cols = [
        "n_pow", "n",
        "build_tree_ms",
        "interp_avg_ms",
        "interp_min_ms",
        "interp_max_ms",
        "interp_ms_per_point",
    ]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["n"]).copy()
    if df.empty:
        print("Error: no valid rows after cleaning CSV.")
        return

    df = df.sort_values(["prime", "n"]) if "prime" in df.columns else df.sort_values(["n"])

    def _plot(metric: str, ylabel: str = None) -> None:
        if metric not in df.columns:
            print(f"Skip: column '{metric}' not found in CSV.")
            return

        plt.figure(figsize=(10, 6))
        if "prime" in df.columns:
            for prime, g in df.groupby("prime"):
                plt.plot(g["n"], g[metric], marker="o", label=f"Prime: {prime}")
            plt.legend()
        else:
            plt.plot(df["n"], df[metric], marker="o")

        plt.xscale("log", base=2)
        if logy:
            plt.yscale("log")

        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.xlabel("n (number of points)")
        plt.ylabel(ylabel or metric)
        plt.title(f"{metric} vs n")
        plt.tight_layout()

        filename = f"{metric}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path)
        print(f"Saved: {save_path}")
        plt.close()

    _plot("build_tree_ms", "build tree time (ms)")
    _plot("interp_avg_ms", "interpolation time avg (ms)")
    _plot("interp_ms_per_point", "interpolation time per point (ms)")


def plot_flint_csv(csv_path: str, logy: bool = True, output_dir: str = "plots") -> None:
    """
    Visualize FLINT bench CSV (lib,mod,k,n,avg_us).
    """
    os.makedirs(output_dir, exist_ok=True)

    df = _load_csv(csv_path)
    if df is None:
        return

    # Remove possible duplicate header rows
    if "mod" in df.columns:
        df = df[df["mod"] != "mod"].copy()

    # Normalize numeric columns
    for c in ["k", "n", "avg_us"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["n", "avg_us"]).copy()
    if df.empty:
        print("Error: no valid rows after cleaning CSV.")
        return

    # Convert to milliseconds for consistency with interp plots
    df["avg_ms"] = df["avg_us"] / 1000.0
    df["ms_per_point"] = df["avg_ms"] / df["n"]

    df = df.sort_values(["mod", "n"]) if "mod" in df.columns else df.sort_values(["n"])

    def _plot(metric: str, ylabel: str) -> None:
        if metric not in df.columns:
            print(f"Skip: column '{metric}' not found in CSV.")
            return

        plt.figure(figsize=(10, 6))
        if "mod" in df.columns:
            for mod, g in df.groupby("mod"):
                plt.plot(g["n"], g[metric], marker="o", label=f"mod: {mod}")
            plt.legend()
        else:
            plt.plot(df["n"], df[metric], marker="o")

        plt.xscale("log", base=2)
        if logy:
            plt.yscale("log")

        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.xlabel("n (number of points)")
        plt.ylabel(ylabel)
        plt.title(f"{metric} vs n")
        plt.tight_layout()

        filename = f"{metric}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path)
        print(f"Saved: {save_path}")
        plt.close()

    _plot("avg_ms", "interpolation time avg (ms)")
    _plot("ms_per_point", "interpolation time per point (ms)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot interpolation benchmark CSV and save figures."
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default="data/bench.csv",
        help="Path to bench CSV (default: data/bench.csv)",
    )
    parser.add_argument(
        "outdir",
        nargs="?",
        default="plots",
        help="Output directory for PNG figures (default: plots/)",
    )
    parser.add_argument(
        "--format",
        choices=["auto", "interp", "flint"],
        default="auto",
        help="Input CSV format (auto-detect, or force interp/flint).",
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="Use log scale on Y axis (default: off unless --logy)",
    )
    parser.add_argument(
        "--linear-y",
        action="store_true",
        help="Force linear Y axis (overrides --logy)",
    )

    args = parser.parse_args()
    logy = args.logy and not args.linear_y
    fmt = args.format

    if fmt == "auto":
        # Heuristic: if build_tree_ms present -> interp, elif avg_us present -> flint, else interp
        df = _load_csv(args.csv)
        if df is None:
            return
        if "build_tree_ms" in df.columns or "interp_avg_ms" in df.columns:
            fmt = "interp"
        elif "avg_us" in df.columns:
            fmt = "flint"
        else:
            fmt = "interp"
    if fmt == "interp":
        plot_interp_csv(args.csv, logy=logy, output_dir=args.outdir)
    else:
        plot_flint_csv(args.csv, logy=logy, output_dir=args.outdir)


if __name__ == "__main__":
    main()
