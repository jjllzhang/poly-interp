import argparse
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")  # no-GUI backend
import matplotlib.pyplot as plt


def plot_interp_csv(csv_path: str, logy: bool = True, output_dir: str = "plots") -> None:
    """
    Visualize benchmark CSV and save figures to disk.
    """
    os.makedirs(output_dir, exist_ok=True)

    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: File '{csv_path}' not found.")
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
    plot_interp_csv(args.csv, logy=logy, output_dir=args.outdir)


if __name__ == "__main__":
    main()
