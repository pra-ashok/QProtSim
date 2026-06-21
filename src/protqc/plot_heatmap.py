"""
Plot a residue-residue contact/weight matrix as a heatmap.

Usage:
    python plot_contact_heatmap.py input.csv output.png [--title "My Title"] [--label "Contact Score"]

Expects a CSV with residue labels as both the index (first column) and the header row,
e.g.:
        ,1-ASN,2-LEU,3-TYR,...
        1-ASN,0.0,0.74,0.30,...
        2-LEU,0.74,0.0,0.74,...
        ...
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_contact_heatmap(
    csv_path: str,
    output_path: str,
    title: str = "Inter-Residue Contact Map",
    cbar_label: str = "Contact Score",
    fixed_vmax: float | None = None,
    font_size: int = 14,
    fig_size: tuple[float, float] = (20, 18),
) -> None:
    df = pd.read_csv(csv_path, index_col=0)

    HEATMAP_CMAP = sns.color_palette("light:darkblue", as_cmap=True)

    # Mask exact zeros so "no contact" cells render as clean white, not the
    # darkest/lightest end of the colormap.
    df_plot = df.replace(0.0, np.nan)

    # Use a fixed vmax (e.g. 1.0) for normalized scores, or None to scale
    # dynamically to the data's own max (useful for raw probabilities).
    vmax = fixed_vmax if fixed_vmax is not None else df_plot.max().max()

    fig, ax = plt.subplots(figsize=fig_size)

    sns.heatmap(
        df_plot,
        ax=ax,
        cmap=HEATMAP_CMAP,
        vmin=0.0,
        vmax=vmax,
        linewidths=0.5,
        linecolor="#e0e0e0",
        square=True,
        annot=True,
        fmt=".2f",
        annot_kws={"size": font_size, "color": "white", "weight": "bold"},
        cbar_kws={"shrink": 0.75, "label": cbar_label, "pad": 0.02},
        mask=df_plot.isna(),
    )

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=font_size)
    cbar.set_label(cbar_label, fontsize=font_size, labelpad=12)

    ax.set_xlabel("Residue", fontsize=font_size, labelpad=12, fontweight="bold")
    ax.set_ylabel("Residue", fontsize=font_size, labelpad=12, fontweight="bold")
    ax.set_title(title, fontsize=font_size, fontweight="bold", pad=18)

    ax.tick_params(axis="x", labelsize=font_size, rotation=45)
    ax.tick_params(axis="y", labelsize=font_size, rotation=0)

    plt.tight_layout()
    plt.savefig(output_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved heatmap to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a residue contact heatmap.")
    parser.add_argument("csv_path", help="Path to input CSV (residue x residue matrix).")
    parser.add_argument("output_path", help="Path to save the output PNG.")
    parser.add_argument("--title", default="Inter-Residue Contact Map")
    parser.add_argument("--label", default="Contact Score", dest="cbar_label")
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Fixed vmax (e.g. 1.0). Omit to scale dynamically to the data's max value.",
    )
    parser.add_argument("--font-size", type=int, default=14)
    args = parser.parse_args()

    plot_contact_heatmap(
        csv_path=args.csv_path,
        output_path=args.output_path,
        title=args.title,
        cbar_label=args.cbar_label,
        fixed_vmax=args.vmax,
        font_size=args.font_size,
    )