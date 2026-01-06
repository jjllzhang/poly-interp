import argparse
import os
import re

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


def _slugify(text: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_-]+", "_", text.strip())


def _basename_label(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith(".csv"):
        base = base[:-4]
    return base or "ops"


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


def plot_interp_compare(csv_a: str, csv_b: str, logy: bool = True, output_dir: str = "plots",
                        label_a: str = None, label_b: str = None) -> None:
    """
    Compare two interp-like CSVs on the same plots (e.g., project vs FLINT).
    """
    os.makedirs(output_dir, exist_ok=True)

    df_a = _load_csv(csv_a)
    df_b = _load_csv(csv_b)
    if df_a is None or df_b is None:
        return

    def _clean(df: pd.DataFrame) -> pd.DataFrame:
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
        return df.dropna(subset=["n"]).copy()

    df_a = _clean(df_a)
    df_b = _clean(df_b)
    if df_a.empty or df_b.empty:
        print("Error: one of the CSVs has no valid rows after cleaning.")
        return

    la = label_a or _basename_label(csv_a)
    lb = label_b or _basename_label(csv_b)
    df_a["source"] = la
    df_b["source"] = lb

    df = pd.concat([df_a, df_b], ignore_index=True)
    df = df.sort_values(["op" if "op" in df.columns else "prime", "n"])

    def _plot(metric: str, ylabel: str) -> None:
        if metric not in df.columns:
            print(f"Skip: column '{metric}' not found in CSV.")
            return

        plt.figure(figsize=(10, 6))
        if "prime" in df.columns:
            for (source, prime), g in df.groupby(["source", "prime"]):
                plt.plot(g["n"], g[metric], marker="o", label=f"{source}-{prime}")
        else:
            for source, g in df.groupby("source"):
                plt.plot(g["n"], g[metric], marker="o", label=source)

        plt.xscale("log", base=2)
        if logy:
            plt.yscale("log")

        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.xlabel("n (number of points)")
        plt.ylabel(ylabel)
        plt.title(f"{metric} vs n (compare)")
        plt.legend()
        plt.tight_layout()

        filename = f"compare_interp_{metric}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path)
        print(f"Saved: {save_path}")
        plt.close()

    _plot("build_tree_ms", "build tree time (ms)")
    _plot("interp_avg_ms", "interpolation time avg (ms)")
    _plot("interp_ms_per_point", "interpolation time per point (ms)")


def plot_ops_csv(csv_path: str, logy: bool = True, output_dir: str = "plots") -> None:
    """
    Visualize ops_bench CSV (prime,op,n_pow,n,avg_ms,min_ms,max_ms,ms_per_item).
    Produces two plots per op: avg_ms vs n, ms_per_item vs n.
    """
    os.makedirs(output_dir, exist_ok=True)

    df = _load_csv(csv_path)
    if df is None:
        return

    # Drop repeated header rows if any
    if "prime" in df.columns:
        df = df[df["prime"] != "prime"].copy()

    num_cols = ["n_pow", "n", "avg_ms", "min_ms", "max_ms", "ms_per_item"]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["n", "op"]).copy()
    if df.empty:
        print("Error: no valid rows after cleaning CSV.")
        return

    df = df.sort_values(["op", "prime", "n"])

    def _plot_op(op: str, metric: str, ylabel: str) -> None:
        if metric not in df.columns:
            print(f"Skip: column '{metric}' not found in CSV.")
            return

        sub = df[df["op"] == op]
        if sub.empty:
            return

        plt.figure(figsize=(10, 6))
        for prime, g in sub.groupby("prime"):
            plt.plot(g["n"], g[metric], marker="o", label=f"prime: {prime}")

        plt.xscale("log", base=2)
        if logy:
            plt.yscale("log")

        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.xlabel("n (number of points)")
        plt.ylabel(ylabel)
        plt.title(f"{op} – {metric} vs n")
        plt.legend()
        plt.tight_layout()

        filename = f"{_slugify(op)}_{metric}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path)
        print(f"Saved: {save_path}")
        plt.close()

    for op in df["op"].unique():
        _plot_op(op, "avg_ms", "time avg (ms)")
        _plot_op(op, "ms_per_item", "time per item (ms)")


def plot_ops_compare(csv_a: str, csv_b: str, logy: bool = True, output_dir: str = "plots",
                     label_a: str = None, label_b: str = None) -> None:
    """
    Plot two ops CSVs on the same chart for comparison (e.g., project vs FLINT).
    """
    os.makedirs(output_dir, exist_ok=True)

    df_a = _load_csv(csv_a)
    df_b = _load_csv(csv_b)
    if df_a is None or df_b is None:
        return

    def _clean(df: pd.DataFrame) -> pd.DataFrame:
        if "prime" in df.columns:
            df = df[df["prime"] != "prime"].copy()
        for c in ["n_pow", "n", "avg_ms", "min_ms", "max_ms", "ms_per_item"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")
        return df.dropna(subset=["n", "op"]).copy()

    df_a = _clean(df_a)
    df_b = _clean(df_b)
    if df_a.empty or df_b.empty:
        print("Error: one of the CSVs has no valid rows after cleaning.")
        return

    la = label_a or _basename_label(csv_a)
    lb = label_b or _basename_label(csv_b)
    df_a["source"] = la
    df_b["source"] = lb

    df = pd.concat([df_a, df_b], ignore_index=True)
    df = df.sort_values(["op", "prime", "source", "n"])

    def _plot_op(op: str, metric: str, ylabel: str) -> None:
        if metric not in df.columns:
            print(f"Skip: column '{metric}' not found in CSV.")
            return
        sub = df[df["op"] == op]
        if sub.empty:
            return

        plt.figure(figsize=(10, 6))
        for (source, prime), g in sub.groupby(["source", "prime"]):
            plt.plot(g["n"], g[metric], marker="o", label=f"{source}-{prime}")

        plt.xscale("log", base=2)
        if logy:
            plt.yscale("log")

        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.xlabel("n (number of points)")
        plt.ylabel(ylabel)
        plt.title(f"{op} – {metric} vs n (compare)")
        plt.legend()
        plt.tight_layout()

        filename = f"compare_{_slugify(op)}_{metric}.png"
        save_path = os.path.join(output_dir, filename)
        plt.savefig(save_path)
        print(f"Saved: {save_path}")
        plt.close()

    for op in df["op"].unique():
        _plot_op(op, "avg_ms", "time avg (ms)")
        _plot_op(op, "ms_per_item", "time per item (ms)")


def plot_flint_csv(csv_path: str, logy: bool = True, output_dir: str = "plots") -> None:
    """
    Visualize FLINT bench CSV (lib,mod,k,n,avg_us).
    """
    os.makedirs(output_dir, exist_ok=True)

    df = _load_csv(csv_path)
    if df is None:
        return

    # If already in interp-like format (aligned with interp_bench), reuse that plotter for consistency.
    if ("build_tree_ms" in df.columns) or ("interp_avg_ms" in df.columns):
        plot_interp_csv(csv_path, logy=logy, output_dir=output_dir)
        return

    # If header is missing (e.g., file appended without header), assign one.
    if "avg_us" not in df.columns:
        if df.shape[1] >= 5:
            df = df.rename(columns={
                df.columns[0]: "lib",
                df.columns[1]: "mod",
                df.columns[2]: "k",
                df.columns[3]: "n",
                df.columns[4]: "avg_us",
            })
    if "n" not in df.columns or "avg_us" not in df.columns:
        print("Error: FLINT CSV missing expected columns 'n'/'avg_us'.")
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
        choices=["auto", "interp", "flint", "ops"],
        default="auto",
        help="Input CSV format (auto-detect, or force interp/flint/ops).",
    )
    parser.add_argument(
        "--compare-ops",
        default=None,
        help="Second ops CSV for comparison (only used with --format=ops or auto-detected ops).",
    )
    parser.add_argument(
        "--compare-label",
        default=None,
        help="Label for the comparison CSV (optional; default uses basename).",
    )
    parser.add_argument(
        "--compare-interp",
        default=None,
        help="Second interp CSV for comparison (only used with --format=interp or auto-detected interp).",
    )
    parser.add_argument(
        "--compare-interp-label",
        default=None,
        help="Label for the comparison interp CSV (optional; default uses basename).",
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
        if "op" in df.columns and "avg_ms" in df.columns:
            fmt = "ops"
        elif "build_tree_ms" in df.columns or "interp_avg_ms" in df.columns:
            fmt = "interp"
        elif "avg_us" in df.columns:
            fmt = "flint"
        else:
            fmt = "interp"
    if fmt == "interp":
        if args.compare_interp:
            plot_interp_compare(
                args.csv,
                args.compare_interp,
                logy=logy,
                output_dir=args.outdir,
                label_a=_basename_label(args.csv),
                label_b=args.compare_interp_label or _basename_label(args.compare_interp),
            )
        else:
            plot_interp_csv(args.csv, logy=logy, output_dir=args.outdir)
    elif fmt == "ops":
        if args.compare_ops:
            plot_ops_compare(args.csv, args.compare_ops, logy=logy,
                             output_dir=args.outdir,
                             label_a=_basename_label(args.csv),
                             label_b=args.compare_label or _basename_label(args.compare_ops))
        else:
            plot_ops_csv(args.csv, logy=logy, output_dir=args.outdir)
    else:
        plot_flint_csv(args.csv, logy=logy, output_dir=args.outdir)


if __name__ == "__main__":
    main()
