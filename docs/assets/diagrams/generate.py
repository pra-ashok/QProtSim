"""Regenerate the documentation's explanatory scientific figures.

Run from any directory with NumPy and Matplotlib installed.
SVG assets are written beside this script. PNG previews are optional via
``--preview-dir`` and support visual verification without raster files in git.
"""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


INK = "#243746"
BLUE = "#28628a"
TEAL = "#267f75"
ORANGE = "#bb622e"
MUTED = "#657583"
PALE = "#f1f5f7"
GRID = "#dce3e7"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "mathtext.fontset": "stix",
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.titleweight": "bold",
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "svg.fonttype": "none",
        "svg.hashsalt": "qprotsim-documentation",
        "savefig.facecolor": "white",
    }
)


def save_figure(fig, name, preview_dir):
    directory = Path(__file__).resolve().parent
    fig.savefig(directory / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    if preview_dir is not None:
        preview_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview_dir / f"{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def box(ax, x, y, width, height, text, color=BLUE, fontsize=11):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.01,rounding_size=0.05",
        linewidth=1.1,
        edgecolor=color,
        facecolor="white",
        zorder=3,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        linespacing=1.6,
        zorder=4,
    )


def arrow(ax, start, end, color=MUTED):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.25,
            color=color,
            shrinkA=5,
            shrinkB=5,
        )
    )


def workflow(preview_dir):
    fig, ax = plt.subplots(figsize=(12.5, 3.7))
    ax.set(xlim=(-0.1, 12.1), ylim=(-0.2, 3.2))
    ax.axis("off")

    box(ax, 0.05, 0.95, 1.75, 1.15, "PDB structure\natomic coordinates")
    box(ax, 2.6, 1.75, 2.4, 0.95, "Relative exposure\n" + r"$\alpha_i=A_i/A_{\max}(a_i)$", TEAL)
    box(ax, 2.6, 0.25, 2.4, 0.95, "Structural contacts\n" + r"$K_{ij}\propto S_{ij}\Phi_{ij}\Lambda_{ij}$", ORANGE)
    box(ax, 5.8, 0.95, 2.5, 1.15, "Quantum circuit\n" + r"$R_Y(\theta_i)$ + ordered pair gates")
    box(ax, 9.1, 0.95, 2.7, 1.15, "Sampled ensemble\n" + r"$Z\in\{+1,-1\}^{S\times N}$")

    arrow(ax, (1.8, 1.63), (2.6, 2.17), TEAL)
    arrow(ax, (1.8, 1.37), (2.6, 0.73), ORANGE)
    arrow(ax, (5.0, 2.17), (5.8, 1.66), TEAL)
    arrow(ax, (5.0, 0.73), (5.8, 1.37), ORANGE)
    arrow(ax, (8.3, 1.53), (9.1, 1.53))

    ax.text(3.8, 2.95, "Local thermodynamic input", ha="center", color=TEAL, fontsize=10)
    ax.text(3.8, -0.05, "Residue coupling graph", ha="center", color=ORANGE, fontsize=10)
    ax.text(7.05, 0.45, r"$\theta_i=2\arcsin\sqrt{p_i}$", ha="center", fontsize=12)
    ax.text(10.45, 0.45, "Energies · ground probabilities\nKL contributions · MI · DEI", ha="center", fontsize=10, linespacing=1.6)
    fig.subplots_adjust(left=0.015, right=0.985, bottom=0.08, top=0.98)
    save_figure(fig, "workflow", preview_dir)


def excitation(preview_dir):
    alpha = np.linspace(0.0, 0.999, 500)
    residues = [
        ("ALA", 1.76 * (1.0 - alpha), BLUE),
        ("GLY", 0.10 * (1.0 - np.maximum(alpha, 0.1)), TEAL),
        ("LYS", 5.61 * alpha, ORANGE),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.7))
    for label, gap, color in residues:
        axes[0].plot(alpha, gap, color=color, label=label, lw=2.4)
        p_exc = 1.0 / (1.0 + np.exp(gap / (8.314462618e-3 * 298.0)))
        axes[1].plot(alpha, p_exc, color=color, label=label, lw=2.4)

    axes[0].set_title("a   Exposure-dependent effective gap", loc="left", fontsize=12, pad=17)
    axes[0].set_ylabel(r"Effective gap $g_i$  (kJ mol$^{-1}$)")
    axes[0].set_ylim(0, 6)
    axes[1].set_title("b   Initial excitation probability", loc="left", fontsize=12, pad=17)
    axes[1].set_ylabel(r"$p_i = [1 + \exp(g_i/RT)]^{-1}$")
    axes[1].set_ylim(0, 0.55)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"Relative solvent exposure $\alpha_i$")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color=GRID, lw=0.7)
        ax.set_axisbelow(True)

    axes[0].annotate(
        "Hydrophilic branch",
        xy=(0.64, 5.61 * 0.64),
        xytext=(0.07, 4.7),
        arrowprops={"arrowstyle": "-", "color": ORANGE},
        color=ORANGE,
        fontsize=10,
    )
    axes[0].annotate(
        "Hydrophobic branch",
        xy=(0.25, 1.76 * 0.75),
        xytext=(0.35, 2.25),
        arrowprops={"arrowstyle": "-", "color": BLUE},
        color=BLUE,
        fontsize=10,
    )
    axes[1].text(0.05, 0.065, r"$T=298\ \mathrm{K}$", transform=axes[1].transAxes, fontsize=12)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.51, -0.005))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.23, top=0.85, wspace=0.28)
    save_figure(fig, "excitation", preview_dir)


def controlled_ry(ax, x, control_y, target_y, label):
    ax.plot([x, x], [control_y, target_y], color=BLUE, lw=1.4, zorder=2)
    ax.add_patch(Circle((x, control_y), 0.045, color=BLUE, zorder=4))
    box(ax, x - 0.28, target_y - 0.17, 0.56, 0.34, label, fontsize=11)


def circuit_information(preview_dir):
    fig = plt.figure(figsize=(12.0, 6.9))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.55, wspace=0.37)
    circuit_ax = fig.add_subplot(grid[0, :])
    circuit_ax.set(xlim=(-0.2, 8.9), ylim=(-0.15, 2.8))
    circuit_ax.axis("off")
    circuit_ax.set_title("a   Local initialization and one strong scheduled pair", loc="left", fontsize=12, pad=8)
    for index, y in enumerate([2.25, 1.35, 0.45]):
        circuit_ax.text(-0.1, y, rf"$i={index}\quad |0\rangle$", ha="left", va="center", fontsize=12)
        circuit_ax.plot([0.9, 8.1], [y, y], color=MUTED, lw=1.1, zorder=0)
        box(circuit_ax, 1.15, y - 0.17, 0.95, 0.34, rf"$R_Y(\theta_{index})$", fontsize=11)
        box(circuit_ax, 7.3, y - 0.17, 0.65, 0.34, r"$Z$", fontsize=11)
        circuit_ax.text(8.22, y, rf"$Z_{{t{index}}}$", ha="left", va="center", fontsize=12)

    controlled_ry(circuit_ax, 3.2, 2.25, 1.35, r"$R_Y(\phi)$")
    controlled_ry(circuit_ax, 4.45, 1.35, 2.25, r"$R_Y(-\phi)$")
    circuit_ax.add_patch(Rectangle((5.15, 0.2), 1.65, 2.3, facecolor=PALE, edgecolor=GRID, lw=1, linestyle="--", zorder=1))
    circuit_ax.text(5.975, 1.35, "Remaining\nscheduled gates", ha="center", va="center", fontsize=10, linespacing=1.6)
    circuit_ax.text(3.85, -0.1, "Source i = 0 → target j = 1", ha="center", fontsize=10, color=BLUE)

    table_ax = fig.add_subplot(grid[1, 0])
    table_ax.axis("off")
    table_ax.set_title("b   Example joint probability distribution", loc="left", fontsize=12, pad=18)
    joint = np.array([[0.70, 0.10], [0.05, 0.15]])
    table = table_ax.table(
        cellText=[["0.70", "0.10"], ["0.05", "0.15"]],
        rowLabels=[r"$R_i=0$", r"$R_i=1$"],
        colLabels=[r"$R_j=0$", r"$R_j=1$"],
        loc="upper center",
        cellLoc="center",
        rowLoc="center",
        bbox=[0.15, 0.38, 0.80, 0.56],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.set_linewidth(0.9)
        cell.set_facecolor(PALE if row == 0 or col == -1 else "white")
    table_ax.text(0.1, 0.20, r"$P(R_i=1)=0.20\qquad P(R_j=1)=0.25$", transform=table_ax.transAxes, fontsize=12)
    table_ax.text(0.1, 0.025, r"$P(R_i=1\mid R_j=1)=0.60$", transform=table_ax.transAxes, fontsize=12, color=TEAL)

    info_ax = fig.add_subplot(grid[1, 1])
    info_ax.axis("off")
    info_ax.set_title("c   Information from the same distribution", loc="left", fontsize=12, pad=18)
    p_i = joint.sum(axis=1)
    p_j = joint.sum(axis=0)
    mi = np.sum(joint * np.log2(joint / np.outer(p_i, p_j)))
    q_i = joint[:, 1] / p_j[1]
    q_j = joint[1, :] / p_i[1]
    dei_i_j = np.sum(q_i * np.log2(q_i / p_i))
    dei_j_i = np.sum(q_j * np.log2(q_j / p_j))
    entries = [
        (0.85, "Mutual information", rf"$I(R_i;R_j)=I(R_j;R_i)={mi:.3f}\ \mathrm{{bits}}$", BLUE),
        (0.49, "Source j → target i", rf"$D_{{ij}}={dei_i_j:.3f}\ \mathrm{{bits}}$", TEAL),
        (0.16, "Source i → target j", rf"$D_{{ji}}={dei_j_i:.3f}\ \mathrm{{bits}}$", ORANGE),
    ]
    for y, title, equation, color in entries:
        info_ax.text(0.0, y, title, transform=info_ax.transAxes, fontsize=10, color=color)
        info_ax.text(0.0, y - 0.15, equation, transform=info_ax.transAxes, fontsize=14)

    fig.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.08)
    save_figure(fig, "circuit-information", preview_dir)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--preview-dir", type=Path)
    arguments = parser.parse_args()
    workflow(arguments.preview_dir)
    excitation(arguments.preview_dir)
    circuit_information(arguments.preview_dir)
