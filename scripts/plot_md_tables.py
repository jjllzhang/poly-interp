#!/usr/bin/env python3
"""
Render tables from a markdown file as PNG images (optionally filtered).

Examples:

    # Render all tables in a markdown file
    python3 scripts/plot_md_tables.py --input data/interp_vs_flint.md --output-dir plots/interp_tables

    # Render only rows where n_pow=20 from ops_vs_flint.md
    python3 scripts/plot_md_tables.py --input data/ops_vs_flint.md --output-dir plots/ops_tables --filter n_pow=20
"""

import argparse
from pathlib import Path
from typing import List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def extract_tables(md_path: Path) -> List[List[str]]:
    """Return a list of raw markdown table blocks (each as list of lines)."""
    lines = md_path.read_text(encoding="utf-8").splitlines()
    tables: List[List[str]] = []
    current: List[str] = []
    for line in lines:
        if line.strip().startswith("|"):
            current.append(line)
        else:
            if current:
                tables.append(current)
                current = []
    if current:
        tables.append(current)
    return tables


def _maybe_number(value: str):
    """Convert numeric strings to numbers; keep '-' or text as-is."""
    val = value.strip()
    if val in {"", "-"}:
        return val
    try:
        num = float(val)
    except ValueError:
        return val
    if np.isclose(num, round(num)):
        return int(round(num))
    return num


def parse_md_table(table_lines: List[str]) -> pd.DataFrame:
    """Parse a markdown table block into a DataFrame."""
    if len(table_lines) < 2:
        raise ValueError("Table does not contain enough rows.")
    header = [h.strip() for h in table_lines[0].strip().strip("|").split("|")]
    data_lines = table_lines[2:]  # Skip separator row
    rows = []
    for line in data_lines:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if not cells or all(c == "" for c in cells):
            continue
        rows.append([_maybe_number(c) for c in cells])
    return pd.DataFrame(rows, columns=header)


def apply_filters(df: pd.DataFrame, filters: Sequence[str]) -> pd.DataFrame:
    """Filter a DataFrame using simple equality expressions like column=value."""
    result = df
    for expr in filters:
        if "=" not in expr:
            print(f"Warning: ignoring filter '{expr}' (expected column=value).")
            continue
        col, raw_val = expr.split("=", 1)
        col = col.strip()
        val = _maybe_number(raw_val)
        if col not in result.columns:
            print(f"Warning: column '{col}' not in table, skipping filter.")
            continue
        result = result[result[col] == val]
    return result


def format_df_for_display(df: pd.DataFrame) -> pd.DataFrame:
    """Format mixed numeric/text data into strings for table rendering."""
    def _fmt(value) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, str):
            return value
        if isinstance(value, (int, np.integer)):
            return f"{int(value)}"
        if isinstance(value, (float, np.floating)):
            if np.isclose(value, round(value)):
                return f"{int(round(value))}"
            return f"{value:.3f}".rstrip("0").rstrip(".")
        return str(value)

    return df.map(_fmt)


def render_table_image(df: pd.DataFrame, title: str, output_path: Path, dpi: int = 300) -> None:
    """Render a DataFrame as a PNG image."""
    display_df = format_df_for_display(df)
    nrows, ncols = display_df.shape
    fig_width = max(10.0, 1.05 * ncols)
    fig_height = max(4.5, 0.55 * nrows + 1.0)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=16)

    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.1, 1.3)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#f0f0f0")
            cell.set_fontsize(10)
            cell.set_text_props(weight="bold")

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render markdown tables as PNG images.")
    parser.add_argument("--input", default="data/interp_vs_flint.md", help="Path to the markdown file.")
    parser.add_argument(
        "--output-dir", default="plots/interp_tables", help="Directory where PNGs will be written."
    )
    parser.add_argument(
        "--table-index",
        type=int,
        default=None,
        help="Optional 1-based table index to render (default: all tables).",
    )
    parser.add_argument(
        "--filter",
        action="append",
        default=[],
        help="Filter expression column=value (can be repeated). Applied after parsing each table.",
    )
    parser.add_argument(
        "--titles",
        nargs="*",
        default=[],
        help="Optional custom titles per table (space-separated). If fewer than tables, remaining use defaults.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output image DPI.")
    args = parser.parse_args()

    md_path = Path(args.input)
    if not md_path.exists():
        raise SystemExit(f"Input file not found: {md_path}")

    tables = extract_tables(md_path)
    if not tables:
        raise SystemExit("No markdown tables found.")

    titles = args.titles if args.titles else []

    for idx, table_lines in enumerate(tables):
        if args.table_index is not None and (idx + 1) != args.table_index:
            continue
        df = parse_md_table(table_lines)
        if args.filter:
            df = apply_filters(df, args.filter)
        if df.empty:
            print(f"Warning: table {idx + 1} has no rows after filtering, skipping.")
            continue

        title = titles[idx] if idx < len(titles) else f"{md_path.stem} – table {idx + 1}"
        outfile = Path(args.output_dir) / f"table_{idx + 1}.png"
        render_table_image(df, title, outfile, dpi=args.dpi)
        print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
