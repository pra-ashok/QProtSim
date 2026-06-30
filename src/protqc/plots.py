from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

HEATMAP_CMAP = sns.color_palette("light:darkblue", as_cmap=True)


def _overlap_area(box_a: tuple[float, float, float, float], box_b: tuple[float, float, float, float]) -> float:
    x0 = max(box_a[0], box_b[0])
    y0 = max(box_a[1], box_b[1])
    x1 = min(box_a[2], box_b[2])
    y1 = min(box_a[3], box_b[3])
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)


def _scatter_label_candidates(base_offset: float = 7.0, rings: int = 4) -> list[tuple[float, float]]:
    candidates: list[tuple[float, float]] = []
    angles = np.deg2rad([45, 135, 315, 225, 0, 90, 180, 270])
    for ring in range(rings):
        radius = base_offset + ring * 6.0
        for angle in angles:
            candidates.append((float(radius * np.cos(angle)), float(radius * np.sin(angle))))
    return candidates


def _annotation_box(
    anchor_px: np.ndarray,
    offset_pts: tuple[float, float],
    label: str,
    fontsize: float,
    dpi: float,
) -> tuple[float, float, float, float]:
    scale = dpi / 72.0
    dx_px = offset_pts[0] * scale
    dy_px = offset_pts[1] * scale
    width = max(12.0, len(str(label)) * fontsize * 0.62 * scale)
    height = fontsize * 1.35 * scale
    pad = 2.0 * scale

    x_anchor = float(anchor_px[0] + dx_px)
    y_anchor = float(anchor_px[1] + dy_px)

    if offset_pts[0] < 0:
        x0, x1 = x_anchor - width - pad, x_anchor + pad
    else:
        x0, x1 = x_anchor - pad, x_anchor + width + pad

    if offset_pts[1] < 0:
        y0, y1 = y_anchor - height - pad, y_anchor + pad
    else:
        y0, y1 = y_anchor - pad, y_anchor + height + pad

    return x0, y0, x1, y1


def _annotate_scatter_labels(
    ax: plt.Axes,
    x_data: Sequence[float],
    y_data: Sequence[float],
    labels: Sequence[str],
    *,
    fontsize: float = 8,
    point_size: float = 70,
    color: str = "black",
    zorder: int = 4,
) -> None:
    """Place scatter labels with deterministic jittered offsets to reduce overlap."""
    fig = ax.figure
    fig.canvas.draw()

    xy = np.column_stack([np.asarray(x_data, dtype=float), np.asarray(y_data, dtype=float)])
    finite_mask = np.isfinite(xy).all(axis=1)
    anchors_px = ax.transData.transform(xy[finite_mask])
    finite_labels = [label for label, keep in zip(labels, finite_mask) if keep]
    if len(finite_labels) == 0:
        return

    dpi = float(fig.dpi)
    point_radius_px = max(4.0, float(np.sqrt(max(point_size, 1.0) / np.pi)) * dpi / 72.0)
    point_boxes = [
        (
            float(px - point_radius_px),
            float(py - point_radius_px),
            float(px + point_radius_px),
            float(py + point_radius_px),
        )
        for px, py in anchors_px
    ]
    axes_box = ax.get_window_extent()
    axes_bounds = (axes_box.x0, axes_box.y0, axes_box.x1, axes_box.y1)
    candidates = _scatter_label_candidates()
    placed_boxes: list[tuple[float, float, float, float]] = []

    for anchor_px, label, xy_point in zip(anchors_px, finite_labels, xy[finite_mask]):
        best_offset = candidates[0]
        best_box = _annotation_box(anchor_px, best_offset, label, fontsize, dpi)
        best_score = float("inf")

        for offset in candidates:
            box = _annotation_box(anchor_px, offset, label, fontsize, dpi)
            label_overlap = sum(_overlap_area(box, other) for other in placed_boxes)
            point_overlap = sum(_overlap_area(box, point_box) for point_box in point_boxes)
            outside = (
                max(0.0, axes_bounds[0] - box[0])
                + max(0.0, axes_bounds[1] - box[1])
                + max(0.0, box[2] - axes_bounds[2])
                + max(0.0, box[3] - axes_bounds[3])
            )
            distance_penalty = 0.04 * (offset[0] * offset[0] + offset[1] * offset[1])
            score = 8.0 * label_overlap + 4.0 * point_overlap + 30.0 * outside + distance_penalty
            if score < best_score:
                best_score = score
                best_offset = offset
                best_box = box

        placed_boxes.append(best_box)
        ax.annotate(
            str(label),
            (float(xy_point[0]), float(xy_point[1])),
            textcoords="offset points",
            xytext=best_offset,
            ha="right" if best_offset[0] < 0 else "left",
            va="top" if best_offset[1] < 0 else "bottom",
            fontsize=fontsize,
            color=color,
            zorder=zorder,
            clip_on=False,
        )

def plot_directional_excitation_information_clustermap(
    csv_path: str,
    save_png: Optional[str] = None,
    save_pdf: Optional[str] = None,
    figsize: tuple[float, float] = (7, 7),
    cmap: str = HEATMAP_CMAP) -> plt.Figure:
    """
    Load DEI (Directional Excitation Information) data from CSV and create a hierarchical clustermap
    using Manhattan distance.

    Parameters
    ----------
    csv_path : str
        Path to the CSV file containing Directional Excitation Information matrix (index_col=0)
    save_png : Optional[str]
        Path to save PNG file. If None, PNG is not saved
    save_pdf : Optional[str]
        Path to save PDF file. If None, PDF is not saved
    figsize : tuple
        Figure size (width, height)
    cmap : str
        Colormap name (default: "magma")

    Returns
    -------
    fig : matplotlib.figure.Figure
        The clustermap figure object
    """
    set_academic_style()

    # Load MI CSV
    mi_df = pd.read_csv(csv_path, index_col=0)

    # Safety: symmetrize and zero diagonal
    mi = mi_df.values
    # mi = 0.5 * (mi + mi.T)
    np.fill_diagonal(mi, 0.0)
    mi_df = pd.DataFrame(mi, index=mi_df.index, columns=mi_df.columns)

    # Normalize MI to [0,1] for stability
    max_mi = mi_df.values.max()
    mi_norm = mi_df / max_mi if max_mi > 0 else mi_df.copy()

    # Manhattan (cityblock) distance clustermap
    g = sns.clustermap(
        mi_norm,
        metric="cityblock",   # Manhattan distance
        method="average",
        cmap=cmap,
        figsize=figsize,
        xticklabels=True,
        yticklabels=True,
    )

    g.fig.suptitle("DEI-based Clustermap (Manhattan / Cityblock distance)", y=1.02)

    # Save to PNG if requested
    if save_png:
        g.fig.savefig(save_png, format="png", dpi=300, bbox_inches="tight")
        print(f"PNG figure saved to: {save_png}")

    # Save to PDF if requested (or derived from PNG path)
    if save_pdf or save_png:
        pdf_path = save_pdf
        if pdf_path is None and save_png:
            pdf_path = str(Path(save_png).with_suffix(".pdf"))
        if pdf_path:
            g.fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
            print(f"PDF figure saved to: {pdf_path}")

    return g.fig  # type: ignore





def set_academic_style() -> None:
    sns.set_theme(
        context="paper",
        style="whitegrid",
        font="DejaVu Serif",
        font_scale=1.15,
        rc={
            "axes.linewidth": 0.8,
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "grid.linewidth": 0.6,
            "grid.alpha": 0.35,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": True,
            "figure.dpi": 600,
            "savefig.dpi": 600,
        },
    )


def _save_figure_with_pdf(fig: plt.Figure, save_path: str, dpi: int = 600) -> None:
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if Path(save_path).suffix.lower() == ".png":
        pdf_path = str(Path(save_path).with_suffix(".pdf"))
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")


def _residue_labels(info_df: pd.DataFrame, results_df: pd.DataFrame) -> list[str]:
    merged = (
        info_df[["res_index"]]
        .merge(results_df[["res_index", "resname"]], on="res_index", how="left")
        .sort_values("res_index")
    )
    return (merged["res_index"].astype(str) + "-" + merged["resname"]).tolist()


def plot_energy_distribution(
    E_total: np.ndarray,
    molecule_name: str,
    bins: int = 1000,
    figsize: tuple[float, float] = (8, 4),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    fig, ax = plt.plot(figsize=figsize)
    sns.histplot(E_total, bins=bins, stat="count", color="#2C3E50", ax=ax)
    #sns.kdeplot(E_total, color="#E67E22", linewidth=1.5, ax=ax)
    ax.set_xlabel("Energy (kJ mol^-1)")
    ax.set_ylabel("Frequency (count)")
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax

def plot_energy_distribution_histogram(
    E_total: np.ndarray,
    molecule_name: str,
    bins: int = 1000,
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None):
    set_academic_style()
    plt.figure(figsize=(8,4))
    bin_edges = np.linspace(min(E_total)-5.0, 0.0, bins + 1)
    counts, bin_edges = np.histogram(E_total, bins=bin_edges, density=False)
    # print("Counts:",len(counts),", Bin edges", len(bin_edges))
    # print("bin edges",bin_edges)
    # print("Counts:",counts)

    plt.bar(bin_edges[:-1], counts)
    plt.plot(bin_edges[:-1], counts)

    plt.xlabel("Energy (kJ·mol⁻¹)")
    plt.ylabel("Frequency (count)")
    if save_path:
        _save_figure_with_pdf(plt.gcf(), save_path, dpi=600)
    return 


def plot_energy_density(
    E_total: np.ndarray,
    molecule_name: str,
    bins: int = 1000,
    figsize: tuple[float, float] = (8, 6),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    fig, ax = plt.subplots(figsize=figsize)
    e_min = float(np.min(E_total))
    e_max = float(np.max(E_total))
    energy_span = max(e_max - e_min, 1.0)
    margin = max(1.0, 0.02 * energy_span)
    bin_edges = np.linspace(e_min - margin, e_max + margin, bins + 1)
    density, bin_edges = np.histogram(E_total, bins=bin_edges, density=True)
    energy_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    ax.fill_between(energy_centers, density, color="#7F93A8", alpha=0.35)
    ax.plot(energy_centers, density, color="#2C3E50", linewidth=1.6)
    ax.set_xlabel("Energy (kJ·mol⁻¹)")
    ax.set_ylabel("Density")
    ax.set_title(f"{molecule_name} energy density")
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax


def plot_energy_sample_scatter(
    E_total: np.ndarray,
    molecule_name: str,
    bins: int = 1000,
    figsize: tuple[float, float] = (8, 6),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    return plot_energy_density(
        E_total,
        molecule_name=molecule_name,
        bins=bins,
        figsize=figsize,
        save_path=save_path,
        fig_num=fig_num,
    )


def plot_dual_axis_comparison(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    results,
    molecule_name: str = "Trp Cage",
    metric1_col: str = "p_ground",
    metric2_col: str = "KL_exc",
    metric1_label: str = r"Stability $P_{ground}$",
    metric2_label: str = r"Energetic Sensitivity $\mathfrak{I}_{\mathrm{exc}}$",
    color1: str = "#2C3E50",
    color2: str = "#E74C3C",
    figsize: tuple[float, float] = (10, 5),
    save_path: str = "trp_cage_comparison_dual_axis.png",
    dpi: int = 600,
    fig_num: Optional[int] = None,
):
    """
    Create a dual-axis bar plot comparing two metrics across residues.

    Parameters
    ----------
    info_df : pd.DataFrame
        DataFrame containing residue information including 'res_index' and metric2_col
    results_df : pd.DataFrame
        DataFrame containing 'resname' column
    results : dict or pd.DataFrame
        Contains metric1_col data
    molecule_name : str
        Name of the molecule for the title
    metric1_col : str
        Column name for first metric (left axis)
    metric2_col : str
        Column name for second metric (right axis)
    metric1_label : str
        Label for first metric (supports LaTeX)
    metric2_label : str
        Label for second metric (supports LaTeX)
    color1 : str
        Color for first metric bars
    color2 : str
        Color for second metric bars
    figsize : tuple
        Figure size (width, height)
    save_path : str
        Path to save the figure
    dpi : int
        Resolution for saved figure

    Returns
    -------
    fig, (ax1, ax2) : matplotlib figure and axes objects
    """
    set_academic_style()

    # Prepare data
    residues = _residue_labels(info_df, results_df)

    # Handle dict or DataFrame for results
    metric1_data = results[metric1_col] if isinstance(results, dict) else results[metric1_col]

    metric2_data = info_df[metric2_col]

    x = np.arange(len(residues))

    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=figsize)

    # Plot bars on primary axis
    width = 0.4
    ax1.bar(
        x - width / 2,
        metric1_data,
        width=width,
        color=color1,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.85,
    )

    ax1.set_xlabel(r"Amino Acid", fontsize=13, fontweight='bold')
    ax1.set_ylabel(metric1_label, fontsize=12, color=color1, fontweight='bold')
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=11)
    ax1.set_ylim(0, max(metric1_data) * 1.15)

    # Create secondary y-axis
    ax2 = ax1.twinx()
    ax2.bar(
        x + width / 2,
        metric2_data,
        width=width,
        color=color2,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.85,
    )
 
    ax2.set_ylabel(f"{metric2_label} (normalised)",
                   fontsize=12, color=color2, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=color2, labelsize=11)
    ax2.set_ylim(0, max(metric2_data) * 1.15)

    # X-axis configuration
    ax1.set_xticks(x)
    ax1.set_xticklabels(residues, rotation=90, ha="right", fontsize=10)
    # Combined legend
    legend_handles = [
        Patch(facecolor=color1, label=metric1_label),
        Patch(facecolor=color2, label=metric2_label),
    ]
    ax1.legend(
        handles=legend_handles,
        #loc="upper right",
        frameon=True,
        fontsize=10,
        bbox_to_anchor=(0.35, 1.2), 
        loc=2, 
        borderaxespad=0.
    )

    # Remove grid
    ax1.grid(False)

    # Add spines for clean look
    sns.despine(ax=ax1, right=False)

    plt.tight_layout()

    # Save figure
    _save_figure_with_pdf(fig, save_path, dpi=dpi)
    print(f"Figure saved to: {save_path}")

    # plt.show()

    return fig, (ax1, ax2)


def plot_residue_information_gain(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    molecule_name: str,
    figsize: tuple[float, float] = (10, 4),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    labels = _residue_labels(info_df, results_df)
    data = pd.DataFrame(
        {
            "Residue": labels,
            "KL_ground": info_df["KL_ground"].values,
            "KL_exc": info_df["KL_exc"].values,
        }
    )
    data = data.melt(id_vars="Residue", var_name="Metric", value_name="KL")

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=data,
        x="Residue",
        y="KL",
        hue="Metric",
        palette=["#2C3E50", "#E74C3C"],
        ax=ax,
        edgecolor="black",
        linewidth=0.4,
    )
    ax.set_xlabel("Residue")
    ax.set_ylabel("Normalised KL contribution")
    ax.tick_params(axis="x", rotation=90, labelsize=9)
    ax.legend(title=None, fontsize=9)
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax


def plot_kl_divergence_scatter(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    molecule_name: str,
    figsize: tuple[float, float] = (6.5, 6),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    del fig_num
    set_academic_style()
    labels = _residue_labels(info_df, results_df)
    x_data = np.asarray(info_df["KL_ground"], dtype=float)
    y_data = np.asarray(info_df["KL_exc"], dtype=float)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(
        x_data,
        y_data,
        s=70,
        color="#2C3E50",
        edgecolor="black",
        linewidth=0.6,
        alpha=0.9,
        zorder=3,
    )

    max_value = float(max(np.max(x_data), np.max(y_data), 1e-12))
    ax.plot(
        [0, max_value * 1.05],
        [0, max_value * 1.05],
        color="#999999",
        linestyle="--",
        linewidth=1.0,
        zorder=1,
    )
    ax.set_xlim(0, max_value * 1.14)
    ax.set_ylim(0, max_value * 1.14)
    _annotate_scatter_labels(ax, x_data, y_data, labels, fontsize=8, point_size=70)
    ax.set_xlabel(r"Normalised $KL_{ground}$")
    ax.set_ylabel(r"Normalised $KL_{exc}$")
    ax.set_title(f"{molecule_name}: KL ground vs. excited")
    ax.grid(False)
    sns.despine(ax=ax)
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax


def plot_kl_divergence_line_comparison(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    molecule_name: str,
    color1: str = "#2C3E50",
    color2: str = "#E74C3C",
    connector_color: str = "#999999",
    connector_alpha: float = 0.6,
    connector_linewidth: float = 0.9,
    figsize: tuple[float, float] = (10, 5),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    del fig_num
    set_academic_style()
    labels = _residue_labels(info_df, results_df)
    ground_data = np.asarray(info_df["KL_ground"], dtype=float)
    exc_data = np.asarray(info_df["KL_exc"], dtype=float)
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=figsize)
    for xi, y_ground, y_exc in zip(x, ground_data, exc_data):
        ax.plot(
            [xi, xi],
            [y_ground, y_exc],
            color=connector_color,
            alpha=connector_alpha,
            linewidth=connector_linewidth,
            linestyle="--",
            zorder=1,
        )

    ax.plot(
        x,
        ground_data,
        color=color1,
        marker="o",
        markersize=5,
        markeredgecolor="black",
        markeredgewidth=0.5,
        linewidth=1.6,
        alpha=0.9,
        label=r"$KL_{ground}$",
        zorder=3,
    )
    ax.plot(
        x,
        exc_data,
        color=color2,
        marker="s",
        markersize=5,
        markeredgecolor="black",
        markeredgewidth=0.5,
        linewidth=1.6,
        alpha=0.9,
        label=r"$KL_{exc}$",
        zorder=3,
    )

    ymax = float(max(np.max(ground_data), np.max(exc_data), 1e-12))
    ax.set_ylim(0, ymax * 1.15)
    ax.set_xlabel("Amino acid", fontsize=13, fontweight="bold")
    ax.set_ylabel("Normalised KL contribution", fontsize=12, fontweight="bold")
    ax.set_title(f"{molecule_name}: Residue KL divergence", fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, ha="right", fontsize=10)
    ax.legend(title=None, frameon=True, fontsize=10)
    ax.grid(False)
    sns.despine(ax=ax)
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax


def plot_stability_vs_information(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    results: dict,
    molecule_name: str,
    figsize: tuple[float, float] = (10, 4),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    labels = _residue_labels(info_df, results_df)
    data = pd.DataFrame(
        {
            "Residue": labels,
            "P_ground": np.asarray(results["p_ground"]),
            "KL_exc": info_df["KL_exc"].values,
        }
    )
    data = data.melt(id_vars="Residue", var_name="Metric", value_name="Value")

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=data,
        x="Residue",
        y="Value",
        hue="Metric",
        palette=["#2C3E50", "#E74C3C"],
        ax=ax,
        edgecolor="black",
        linewidth=0.4,
    )
    ax.set_xlabel("Amino acid")
    ax.set_ylabel("Ground-state probability / normalised KL contribution")
    ax.tick_params(axis="x", rotation=90, labelsize=9)
    ax.legend(title=None, fontsize=9)
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax


def plot_directional_excitation_information_heatmap(
    DEI: np.ndarray,
    labels: Sequence[str],
    figsize: tuple[float, float] = (6, 5),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        DEI,
        cmap=HEATMAP_CMAP,
        ax=ax,
        cbar_kws={"label": "DEI (bits)"},
        xticklabels=labels,
        yticklabels=labels,
    )
    # ax.xaxis.set_label_position("top")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Residue index")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", top=True, bottom=False, rotation=90, labelsize=7)
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax


def plot_pairwise_matrix_heatmap(
    matrix: np.ndarray,
    labels: Sequence[str],
    title: str,
    colorbar_label: str,
    cmap: str = HEATMAP_CMAP,
    fig_size: Optional[tuple[float, float]] = None,
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
    font_size: int = 10,
    annotate: bool = False,
):
    del fig_num
    set_academic_style()
    matrix = np.asarray(matrix, dtype=float)
    n_res = len(labels)
    if fig_size is None:
        side = min(max(0.32 * max(n_res, 1) + 3.2, 6.0), 12.0)
        fig_size = (side, side)

    if isinstance(cmap, str):
        heatmap_cmap = plt.get_cmap(cmap).copy()
    elif hasattr(cmap, "copy"):
        heatmap_cmap = cmap.copy()
    else:
        heatmap_cmap = cmap
    if hasattr(heatmap_cmap, "set_bad"):
        heatmap_cmap.set_bad(color="white")

    data_max = float(np.nanmax(matrix)) if matrix.size else 0.0
    zero_mask = np.isclose(matrix, 0.0) if data_max > 0 else np.zeros_like(matrix, dtype=bool)
    plot_matrix = np.ma.masked_where(zero_mask, matrix)
    vmax = data_max if data_max > 0 else 1.0

    fig, ax = plt.subplots(figsize=fig_size)
    hm = sns.heatmap(
        plot_matrix,
        cmap=heatmap_cmap,
        ax=ax,
        cbar_kws={"label": colorbar_label, "shrink": 0.78, "pad": 0.02},
        linewidths=0.35,
        linecolor="#f0f0f0",
        square=True,
        annot=annotate,
        fmt=".2f",
        annot_kws={"size": max(font_size - 2, 7), "color": "black", "weight": "normal"},
        xticklabels=labels,
        yticklabels=labels,
        vmin=0.0,
        vmax=vmax,
        mask=zero_mask,
    )
    ax.set_facecolor("white")
    ax.set_title(title, fontsize=font_size + 3, pad=14, fontweight="bold")
    ax.set_xlabel("Residue", fontsize=font_size + 1, labelpad=10)
    ax.set_ylabel("Residue", fontsize=font_size + 1, labelpad=10)

    cbar = hm.collections[0].colorbar
    if cbar is not None:
        cbar.set_label(colorbar_label, fontsize=font_size + 1, labelpad=10)
        cbar.ax.tick_params(labelsize=font_size)

    ax.tick_params(axis="x", top=False, bottom=True, rotation=90, labelsize=font_size)
    ax.tick_params(axis="y", rotation=0, labelsize=font_size)
    ax.set_xticklabels(ax.get_xticklabels(), ha="right", rotation_mode="anchor")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)
        spine.set_color("#bdbdbd")
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax

# def plot_contact_heatmap(
#     csv_path: str,
#     output_path: str,
#     title: str = "Inter-Residue Contact Map",
#     cbar_label: str = "Contact Score",
#     fixed_vmax: float | None = None,
#     font_size: int = 14,
#     fig_size: tuple[float, float] = (20, 18),
# ) -> None:
#     df = pd.read_csv(csv_path, index_col=0)

#     HEATMAP_CMAP = sns.color_palette("light:darkblue", as_cmap=True)

#     # Mask exact zeros so "no contact" cells render as clean white, not the
#     # darkest/lightest end of the colormap.
#     df_plot = df.replace(0.0, np.nan)

#     # Use a fixed vmax (e.g. 1.0) for normalized scores, or None to scale
#     # dynamically to the data's own max (useful for raw probabilities).
#     vmax = fixed_vmax if fixed_vmax is not None else df_plot.max().max()

#     fig, ax = plt.subplots(figsize=fig_size)

#     sns.heatmap(
#         df_plot,
#         ax=ax,
#         cmap=HEATMAP_CMAP,
#         vmin=0.0,
#         vmax=vmax,
#         linewidths=0.5,
#         linecolor="#e0e0e0",
#         square=True,
#         annot=True,
#         fmt=".2f",
#         annot_kws={"size": font_size, "color": "white", "weight": "bold"},
#         cbar_kws={"shrink": 0.75, "label": cbar_label, "pad": 0.02},
#         mask=df_plot.isna(),
#     )

#     cbar = ax.collections[0].colorbar
#     cbar.ax.tick_params(labelsize=font_size)
#     cbar.set_label(cbar_label, fontsize=font_size, labelpad=12)

#     ax.set_xlabel("Residue", fontsize=font_size, labelpad=12, fontweight="bold")
#     ax.set_ylabel("Residue", fontsize=font_size, labelpad=12, fontweight="bold")
#     ax.set_title(title, fontsize=font_size, fontweight="bold", pad=18)

#     ax.tick_params(axis="x", labelsize=font_size, rotation=45)
#     ax.tick_params(axis="y", labelsize=font_size, rotation=0)

#     plt.tight_layout()
#     plt.savefig(output_path, dpi=180, bbox_inches="tight", facecolor="white")
#     plt.close(fig)
#     print(f"Saved heatmap to {output_path}")



def plot_joint_state_frequencies(
    matrices: Sequence[np.ndarray],
    titles: Sequence[str],
    labels: Sequence[str],
    figsize: tuple[float, float] = (10, 8),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    for ax, mat, title in zip(axes.flat, matrices, titles):
        # norm = mat / np.max(mat) if np.max(mat) > 0 else mat
        sns.heatmap(
            mat,
            cmap=HEATMAP_CMAP,
            ax=ax,
            cbar_kws={"label": "Normalised frequency"},
            xticklabels=labels,
            yticklabels=labels,
        )
        ax.set_title(title)
        ax.set_xlabel("Residue index")
        ax.set_ylabel("Residue index")
        ax.tick_params(axis="x", rotation=90, labelsize=7)
        ax.tick_params(axis="y", labelsize=7)

    fig.tight_layout(rect=[0, 0, 1, 1])
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, axes


def plot_ground_state_comparison(
    results_df: pd.DataFrame,
    results: dict,
    molecule_name: str,
    figsize: tuple[float, float] = (10, 4),
    save_path: Optional[str] = None,
    fig_num: Optional[int] = None,
):
    set_academic_style()
    labels = (results_df["res_index"].astype(str) + "-" + results_df["resname"]).tolist()
    after_data = np.asarray(results["p_ground"], dtype=float)
    before_data = 1.0 - np.asarray(results["p_exc"], dtype=float)
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=figsize)
    for xi, before_val, after_val in zip(x, before_data, after_data):
        ax.plot(
            [xi, xi],
            [before_val, after_val],
            color="#999999",
            alpha=0.6,
            linewidth=0.9,
            linestyle="--",
            zorder=1,
        )

    ax.plot(
        x,
        after_data,
        color="#2C3E50",
        marker="o",
        markersize=5,
        markeredgecolor="black",
        markeredgewidth=0.5,
        linewidth=1.6,
        alpha=0.9,
        label="After Entanglement",
        zorder=3,
    )
    ax.plot(
        x,
        before_data,
        color="#E74C3C",
        marker="s",
        markersize=5,
        markeredgecolor="black",
        markeredgewidth=0.5,
        linewidth=1.6,
        alpha=0.9,
        label="Before Entanglement",
        zorder=3,
    )

    ymax = float(max(np.max(before_data), np.max(after_data), 1e-12))
    ax.set_ylim(0, min(1.05, ymax * 1.15))
    ax.set_xlabel("Amino acid")
    ax.set_ylabel("P_ground")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, ha="right", fontsize=9)
    ax.legend(title=None, fontsize=9, frameon=True)
    ax.grid(False)
    sns.despine(ax=ax)
    fig.tight_layout()
    if save_path:
        _save_figure_with_pdf(fig, save_path, dpi=600)
    return fig, ax

"""
Dual-axis line plot comparing two metrics across residues, with a vertical
connector line at each residue linking the metric1 point (left axis) to the
metric2 point (right axis).
"""

def plot_dual_axis_line_comparison(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    results,
    molecule_name: str = "Trp Cage",
    metric1_col: str = "p_ground",
    metric2_col: str = "KL_exc",
    metric1_label: str = r"Stability $P_{ground}$",
    metric2_label: str = r"Energetic Sensitivity $\mathfrak{I}_{\mathrm{exc}}$",
    color1: str = "#2C3E50",
    color2: str = "#E74C3C",
    connector_color: str = "#999999",
    connector_alpha: float = 0.6,
    connector_linewidth: float = 0.9,
    figsize: tuple[float, float] = (10, 5),
    save_path: str = "trp_cage_comparison_dual_axis_line.png",
    dpi: int = 600,
    fig_num: Optional[int] = None,
):
    """
    Create a dual-axis line plot comparing two metrics across residues, with
    a vertical line at each residue connecting the metric1 value (left axis)
    to the metric2 value (right axis).

    Parameters
    ----------
    info_df : pd.DataFrame
        DataFrame containing residue information including 'res_index' and metric2_col
    results_df : pd.DataFrame
        DataFrame containing 'resname' column
    results : dict or pd.DataFrame
        Contains metric1_col data
    molecule_name : str
        Name of the molecule for the title
    metric1_col : str
        Column name for first metric (left axis)
    metric2_col : str
        Column name for second metric (right axis)
    metric1_label : str
        Label for first metric (supports LaTeX)
    metric2_label : str
        Label for second metric (supports LaTeX)
    color1 : str
        Color for first metric line/markers
    color2 : str
        Color for second metric line/markers
    connector_color : str
        Color for the vertical connector lines linking the two metrics at each residue
    connector_alpha : float
        Alpha (transparency) for the connector lines
    connector_linewidth : float
        Linewidth for the connector lines
    figsize : tuple
        Figure size (width, height)
    save_path : str
        Path to save the figure
    dpi : int
        Resolution for saved figure

    Returns
    -------
    fig, (ax1, ax2) : matplotlib figure and axes objects
    """
    set_academic_style()

    # Prepare data
    residues = _residue_labels(info_df, results_df)

    # Handle dict or DataFrame for results
    metric1_data = np.asarray(results[metric1_col])
    metric2_data = np.asarray(info_df[metric2_col])

    x = np.arange(len(residues))

    # Create figure and primary axis
    fig, ax1 = plt.subplots(figsize=figsize)

    # Plot metric1 line on primary axis
    ax1.plot(
        x,
        metric1_data,
        color=color1,
        marker="o",
        markersize=5,
        markeredgecolor="black",
        markeredgewidth=0.5,
        linewidth=1.6,
        alpha=0.9,
        zorder=3,
    )

    ax1.set_xlabel(r"Amino Acid", fontsize=13, fontweight="bold")
    ax1.set_ylabel(metric1_label, fontsize=12, color=color1, fontweight="bold")
    ax1.tick_params(axis="y", labelcolor=color1, labelsize=11)
    ax1.set_ylim(0, max(metric1_data) * 1.15)

    # Create secondary y-axis
    ax2 = ax1.twinx()
    ax2.plot(
        x,
        metric2_data,
        color=color2,
        marker="s",
        markersize=5,
        markeredgecolor="black",
        markeredgewidth=0.5,
        linewidth=1.6,
        alpha=0.9,
        zorder=3,
    )

    ax2.set_ylabel(f"{metric2_label} (normalised)", fontsize=12, color=color2, fontweight="bold")
    ax2.tick_params(axis="y", labelcolor=color2, labelsize=11)
    ax2.set_ylim(0, max(metric2_data) * 1.15)

    # Vertical connector lines: from the metric1 point (ax1 data coords) to
    # the metric2 point (ax2 data coords) at the same residue. Since the two
    # axes have different y-scales, each connector is drawn directly on the
    # axis whose transform spans both — here we draw on ax2 using a
    # transform that mixes ax1's y-data with ax2's y-data via display coords.
    for xi, y1, y2 in zip(x, metric1_data, metric2_data):
        # Convert the ax1 point to display coordinates, then to ax2 data coordinates
        disp_pt = ax1.transData.transform((xi, y1))
        y1_in_ax2_coords = ax2.transData.inverted().transform(disp_pt)[1]

        ax2.plot(
            [xi, xi],
            [y1_in_ax2_coords, y2],
            color=connector_color,
            alpha=connector_alpha,
            linewidth=connector_linewidth,
            linestyle="--",
            zorder=1,
        )

    # X-axis configuration
    ax1.set_xticks(x)
    ax1.set_xticklabels(residues, rotation=90, ha="right", fontsize=10)

    # Combined legend
    legend_handles = [
        Line2D([0], [0], color=color1, marker="o", markeredgecolor="black",
               markeredgewidth=0.5, linewidth=1.6, label=metric1_label),
        Line2D([0], [0], color=color2, marker="s", markeredgecolor="black",
               markeredgewidth=0.5, linewidth=1.6, label=metric2_label),
    ]
    ax1.legend(
        handles=legend_handles,
        frameon=True,
        fontsize=10,
        bbox_to_anchor=(0.35, 1.2),
        loc=2,
        borderaxespad=0.0,
    )

    # Remove grid
    ax1.grid(False)

    # Add spines for clean look
    sns.despine(ax=ax1, right=False)

    plt.tight_layout()

    # Save figure
    _save_figure_with_pdf(fig, save_path, dpi=dpi)
    print(f"Figure saved to: {save_path}")

    return fig, (ax1, ax2)

"""
Quadrant scatter plot: Energetic Sensitivity vs Stability, with dotted
reference lines splitting the plot into four quadrants.
"""
def plot_stability_vs_sensitivity_quadrants(
    info_df: pd.DataFrame,
    results_df: pd.DataFrame,
    results,
    molecule_name: str = "Trp Cage",
    metric1_col: str = "p_ground",
    metric2_col: str = "KL_exc",
    metric1_label: str = r"Stability $P_{ground}$",
    metric2_label: str = r"Energetic Sensitivity $\mathfrak{I}_{\mathrm{exc}}$",
    quadrant_method: str = "mean",
    x_threshold: Optional[float] = None,
    y_threshold: Optional[float] = None,
    point_color: str = "#2C3E50",
    quadrant_line_color: str = "#999999",
    annotate_residues: bool = True,
    figsize: tuple[float, float] = (7, 7),
    save_path: str = "trp_cage_stability_vs_sensitivity_quadrants.png",
    dpi: int = 600,
    fig_num: Optional[int] = None,
):
    """
    Scatter plot of Energetic Sensitivity (y-axis) vs Stability (x-axis) for
    each residue, with dotted threshold lines dividing the plot into four
    quadrants.

    Parameters
    ----------
    info_df : pd.DataFrame
        DataFrame containing residue information including 'res_index' and metric2_col
    results_df : pd.DataFrame
        DataFrame containing 'resname' column
    results : dict or pd.DataFrame
        Contains metric1_col data
    molecule_name : str
        Name of the molecule for the title
    metric1_col : str
        Column name for stability metric (x-axis)
    metric2_col : str
        Column name for sensitivity metric (y-axis)
    metric1_label : str
        Label for x-axis (supports LaTeX)
    metric2_label : str
        Label for y-axis (supports LaTeX)
    quadrant_method : str
        How to position the quadrant lines when explicit thresholds aren't
        given: "mean" or "median" of each metric.
    x_threshold, y_threshold : float, optional
        Explicit values for the vertical/horizontal quadrant lines. If None,
        computed from the data using quadrant_method.
    point_color : str
        Color for the scatter points
    quadrant_line_color : str
        Color for the dotted quadrant divider lines
    annotate_residues : bool
        If True, label each point with its residue (e.g. "5-GLN")
    figsize : tuple
        Figure size (width, height)
    save_path : str
        Path to save the figure
    dpi : int
        Resolution for saved figure

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
    """
    set_academic_style()

    residues = _residue_labels(info_df, results_df)

    x_data = np.asarray(results[metric1_col])
    y_data = np.asarray(info_df[metric2_col])

    if x_threshold is None:
        x_threshold = np.mean(x_data) if quadrant_method == "mean" else np.median(x_data)
    if y_threshold is None:
        y_threshold = np.mean(y_data) if quadrant_method == "mean" else np.median(y_data)

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        x_data,
        y_data,
        s=70,
        color=point_color,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.85,
        zorder=3,
    )

    if annotate_residues:
        x_min = float(min(np.min(x_data), x_threshold))
        x_max = float(max(np.max(x_data), x_threshold))
        y_min = float(min(np.min(y_data), y_threshold))
        y_max = float(max(np.max(y_data), y_threshold))
        x_pad = max((x_max - x_min) * 0.10, 1e-3)
        y_pad = max((y_max - y_min) * 0.10, 1e-3)
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)
        _annotate_scatter_labels(
            ax,
            x_data,
            y_data,
            residues,
            fontsize=8,
            point_size=70,
        )

    # Quadrant divider lines
    ax.axvline(x_threshold, color=quadrant_line_color, linestyle="--", linewidth=1.2, zorder=1)
    ax.axhline(y_threshold, color=quadrant_line_color, linestyle="--", linewidth=1.2, zorder=1)

    ax.set_xlabel(metric1_label, fontsize=13, fontweight="bold")
    ax.set_ylabel(f"{metric2_label} (normalised)", fontsize=13, fontweight="bold")
    ax.set_title(f"{molecule_name}: Stability vs. Energetic Sensitivity", fontsize=13, fontweight="bold")

    ax.grid(False)
    sns.despine(ax=ax)

    plt.tight_layout()

    _save_figure_with_pdf(fig, save_path, dpi=dpi)
    print(f"Figure saved to: {save_path}")

    return fig, ax
