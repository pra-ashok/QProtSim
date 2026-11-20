from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional, ParamSpec
import gc
import matplotlib.pyplot as plt
import psutil
import seaborn as sns
import numpy as np
import pandas as pd
import yaml
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib.utils import ImageReader
from reportlab.platypus import Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from qprotsim.circuit_qk import Config
from qprotsim.info_theory import (
    analyse_motif_allostery,
    directional_excitation_information,
    joint_state_matrices,
    mutual_information_matrix,
    residue_information_gain,
)
from qprotsim.pipeline import export_entanglement_parameters, export_residue_statistics, run_pipeline
from qprotsim.plots import (
    plot_dual_axis_line_comparison,
    plot_energy_density,
    plot_energy_distribution_histogram,
    plot_ground_state_comparison,
    plot_kl_divergence_line_comparison,
    plot_kl_divergence_scatter,
    plot_joint_state_frequencies,
    plot_pairwise_matrix_heatmap,
    plot_directional_excitation_information_clustermap,
    plot_directional_excitation_information_heatmap,
    plot_stability_vs_sensitivity_quadrants,
)
from qprotsim.structure_viz import get_structure_image, render_structure_pymol
from memory_profiler import profile

_P = ParamSpec("_P")

@dataclass
class AnalysisOutputs:
    results: dict
    results_df: pd.DataFrame
    info_df: pd.DataFrame
    dei_matrix: np.ndarray
    mi_matrix: np.ndarray
    joint_counts: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    plot_paths: dict[str, Path]
    fig_numbers: dict[str, int]
    run_stamp: str


def load_yaml_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError("Config YAML must be a mapping of keys to values.")
    return data


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _figure_path(output_dir: Path, run_stamp: str, fig_num: int, name: str, suffix: str = ".png") -> Path:
    del run_stamp
    return output_dir / f"fig{fig_num:02d}_{name}{suffix}"


def _figure_label(fig_num: int) -> str:
    return f"Figure {fig_num:02d}"


def _render_report_figure(
    plot: Callable[_P, object], *args: _P.args, **kwargs: _P.kwargs
) -> None:
    existing_figures = set(plt.get_fignums())
    try:
        plot(*args, **kwargs)
    finally:
        for number in set(plt.get_fignums()) - existing_figures:
            plt.close(number)
        gc.collect()


def _format_run_stamp(run_stamp: str) -> str:
    try:
        return datetime.strptime(run_stamp, "%Y%m%d_%H%M%S").strftime("%Y-%m-%d %H:%M:%S")
    except ValueError:
        return run_stamp


def _matrix_from_pair_values(
    pair_values: dict[tuple[int, int], float],
    n_res: int,
) -> np.ndarray:
    matrix = np.zeros((n_res, n_res), dtype=float)
    for (i, j), value in pair_values.items():
        if 0 <= i < n_res and 0 <= j < n_res and i != j:
            v = float(value)
            matrix[i, j] = v
            matrix[j, i] = v
    np.fill_diagonal(matrix, 0.0)
    return matrix


def _safe_artifact_name(name: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in name).strip("_")


def _resolve_output_paths(
    config_path: Path,
    output_dir: Optional[Path],
    config_output_dir: Optional[Path],
    report_filename: Optional[str],
    config_report_filename: Optional[str],
    molecule_name: str,
    run_stamp: Optional[str] = None,
) -> tuple[Path, Path]:
    """
    Resolve output paths with priority: function param > config file > default
    
    Priority order:
    1. Function parameter (output_dir) - highest priority
    2. Config file value (config_output_dir)
    3. Default location (config_path.parent / "qprotsim_outputs")
    
    Report file is ALWAYS placed at the root of the output directory.
    Report filename priority:
    1. Function parameter (report_filename)
    2. Config file value (config_report_filename)
    3. Default: report.pdf
    """
    # Resolve output directory
    if output_dir is not None:
        out_dir = output_dir
    elif config_output_dir is not None:
        out_dir = config_output_dir
    else:
        out_dir = config_path.parent / "qprotsim_outputs"
    
    if report_filename is not None:
        rpt_filename = report_filename
    elif config_report_filename is not None:
        rpt_filename = config_report_filename
    else:
        rpt_filename = "report.pdf"
    
    # Report is always at the root of output_dir
    rpt_path = out_dir / rpt_filename
    
    return out_dir, rpt_path


def _parse_chain_value(data: dict) -> Optional[str]:
    """
    Parse chain value from YAML config.

    Returns None for null/empty/"--" values to indicate no chain filter.
    """
    raw_chain = data.get("chain", "--")
    if raw_chain is None:
        return None
    chain = str(raw_chain).strip()
    if not chain or chain == "--":
        return None
    return chain


def _config_from_data(data: dict) -> tuple[Config, str, Optional[str], int, Optional[Path], Optional[str]]:
    if "pdb_path" not in data and "pdb_code" not in data:
        raise ValueError("Missing required key: either 'pdb_path' or 'pdb_code' must be provided")
    if "molecule_name" not in data:
        raise ValueError("Missing required key: molecule_name")

    # Handle PDB input - either as file path or PDB code to download
    pdb_path = None
    if "pdb_path" in data:
        # User provided a file path - use it directly
        pdb_path = str(data["pdb_path"])
        if not Path(pdb_path).exists():
            raise FileNotFoundError(
                f"\n{'='*70}\n"
                f"ERROR: PDB file not found\n"
                f"{'='*70}\n"
                f"The specified PDB file does not exist:\n"
                f"  {pdb_path}\n\n"
                f"Please check:\n"
                f"  1. The file path in your YAML config is correct\n"
                f"  2. The file exists at the specified location\n"
                f"  3. You have permission to access the file\n"
                f"{'='*70}\n"
            )
    elif "pdb_code" in data:
        # User provided a PDB code - download if needed
        from qprotsim.helper import fetch_pdb
        pdb_code = str(data["pdb_code"])
        output_dir = data.get("download_dir", None)
        
        print(f"\n[INFO] PDB code provided: {pdb_code}")
        print("[INFO] Checking if download is needed...")
        
        pdb_path = str(fetch_pdb(
            pdb_code=pdb_code,
            output_dir=output_dir,
            force_download=data.get("force_download", False),
            verify_checksum=data.get("verify_checksum", True)
        ))
        
        print(f"[INFO] Using PDB file: {pdb_path}\n")
    
    molecule_name = str(data["molecule_name"])
    chain = _parse_chain_value(data)
    bins = int(data.get("bins", 1000))
    
    # Extract output paths from config if provided
    config_output_dir = Path(data["output_dir"]) if "output_dir" in data else None
    config_report_filename = data.get("report_file_name", None)  # Just the filename, not full path

    cfg = Config(
        pdb_path=pdb_path,
        protein_only=bool(data.get("protein_only", True)),
        shots=int(data.get("shots", 100_000)),
        pl_device=str(data.get("pl_device", "default.qubit")),
    )
    return cfg, molecule_name, chain, bins, config_output_dir, config_report_filename

@profile
def run_analysis(
    cfg: Config,
    molecule_name: str,
    chain: Optional[str],
    output_dir: Path,
    bins: int = 1000,
    run_stamp: Optional[str] = None,
) -> AnalysisOutputs:
    """Sample the model and export tables and figures.

    Returns numerical results, information matrices, and generated artifact
    paths. Call generate_report with the returned object to write a PDF.
    bins controls energy plots. Residue information uses its own default
    histogram resolution. run_stamp supplies report metadata, while artifact
    filenames remain stable within output_dir.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    slug = molecule_name.lower().replace(" ", "_")
    run_stamp = run_stamp or _timestamp()

    plot_paths: dict[str, Path] = {}
    fig_numbers: dict[str, int] = {}
    fig_num = 1

    # Generate structure visualization first to stabilize figure numbering.
    # Try to extract PDB code from path (e.g., "1l2y.pdb" -> "1L2Y")
    pdb_code = None
    pdb_filename = Path(cfg.pdb_path).stem.upper()
    if len(pdb_filename) == 4 and pdb_filename.isalnum():
        pdb_code = pdb_filename

    # Get RCSB structure image (downloaded)
    structure_image = get_structure_image(
        pdb_path=Path(cfg.pdb_path),
        output_dir=output_dir,
        pdb_code=pdb_code,
        chain=chain,
        prefer_download=True,
        output_path=_figure_path(output_dir, run_stamp, fig_num, "structure_rcsb"),
    )
    if structure_image:
        plot_paths["structure_rcsb"] = structure_image
        fig_numbers["structure_rcsb"] = fig_num
        fig_num += 1

    # Render with PyMOL (custom rendering)
    pymol_image = render_structure_pymol(
        pdb_path=Path(cfg.pdb_path),
        output_path=_figure_path(output_dir, run_stamp, fig_num, "structure_pymol"),
        width=800,
        height=600,
        style="cartoon",
        color_scheme="spectrum",
        ray_trace=True,
    )
    if pymol_image:
        plot_paths["structure_pymol"] = pymol_image
        fig_numbers["structure_pymol"] = fig_num
        fig_num += 1

    results = run_pipeline(cfg, molecule_name, chain=chain, plot_energy=False, bins=bins)
    results_df = export_residue_statistics(
        results,
        out_csv=str(output_dir / f"{slug}_residue_stats.csv"),
    )
    export_entanglement_parameters(
        results,
        out_csv=str(output_dir / f"{slug}_entanglement_parameters.csv"),
    )
    results_df.sort_values("res_index", inplace=True)
    results_df.reset_index(drop=True, inplace=True)

    info_df = residue_information_gain(results["Z_samples"], results["E"])
    info_df.to_csv(output_dir / f"{slug}_residue_information_gain.csv", index=False)
    info_df = info_df.merge(results_df[["res_index", "resname"]], on="res_index", how="left")
    info_df.sort_values("res_index", inplace=True)
    info_df.reset_index(drop=True, inplace=True)

    plot_paths["energy_distribution"] = _figure_path(
        output_dir, run_stamp, fig_num, "energy_distribution"
    )
    _render_report_figure(
        plot_energy_distribution_histogram,
        results["E"],
        molecule_name=_figure_label(fig_num),
        bins=bins,
        save_path=str(plot_paths["energy_distribution"]),
        fig_num=fig_num,
    )
    fig_numbers["energy_distribution"] = fig_num
    fig_num += 1

    plot_paths["energy_density"] = _figure_path(output_dir, run_stamp, fig_num, "energy_density")
    _render_report_figure(
        plot_energy_density,
        results["E"],
        molecule_name=_figure_label(fig_num),
        bins=bins,
        save_path=str(plot_paths["energy_density"]),
        fig_num=fig_num,
    )
    fig_numbers["energy_density"] = fig_num
    fig_num += 1

    plot_paths["kl_divergence"] = _figure_path(output_dir, run_stamp, fig_num, "kl_divergence")
    _render_report_figure(
        plot_kl_divergence_scatter,
        info_df,
        results_df,
        molecule_name=_figure_label(fig_num),
        save_path=str(plot_paths["kl_divergence"]),
        fig_num=fig_num,
    )
    fig_numbers["kl_divergence"] = fig_num
    fig_num += 1

    plot_paths["kl_divergence_line"] = _figure_path(
        output_dir, run_stamp, fig_num, "kl_divergence_line"
    )
    _render_report_figure(
        plot_kl_divergence_line_comparison,
        info_df,
        results_df,
        molecule_name=_figure_label(fig_num),
        save_path=str(plot_paths["kl_divergence_line"]),
        fig_num=fig_num,
    )
    fig_numbers["kl_divergence_line"] = fig_num
    fig_num += 1

    plot_paths["stability_vs_sensitivity_quadrants"] = (
        _figure_path(output_dir, run_stamp, fig_num, "stability_vs_sensitivity_quadrants")
    )
    plot_paths["stability_vs_info"] = plot_paths["stability_vs_sensitivity_quadrants"]
    _render_report_figure(
        plot_stability_vs_sensitivity_quadrants,
        info_df,
        results_df,
        results,
        molecule_name=_figure_label(fig_num),
        save_path=str(plot_paths["stability_vs_sensitivity_quadrants"]),
        fig_num=fig_num,
    )
    fig_numbers["stability_vs_sensitivity_quadrants"] = fig_num
    fig_numbers["stability_vs_info"] = fig_num
    fig_num += 1

    plot_paths["dual_axis"] = _figure_path(output_dir, run_stamp, fig_num, "dual_axis_comparison")
    _render_report_figure(
        plot_dual_axis_line_comparison,
        info_df=info_df,
        results_df=results_df,
        results=results,
        molecule_name=_figure_label(fig_num),
        save_path=str(plot_paths["dual_axis"]),
        fig_num=fig_num,
    )
    fig_numbers["dual_axis"] = fig_num
    fig_num += 1

    labels = [f"{i}-{n}" for i, n in zip(results_df["res_index"], results_df["resname"])]
    z_samples_path = output_dir / f"{slug}_z_samples.csv"
    z_samples_df = pd.DataFrame(results["Z_samples"], columns=labels)
    z_samples_df.index.name = "shot"
    z_samples_df.to_csv(z_samples_path)

    n_res = len(labels)
    kij_matrix = _matrix_from_pair_values(results.get("kijs", {}), n_res)
    pd.DataFrame(kij_matrix, index=labels, columns=labels).to_csv(
        output_dir / f"{slug}_kij_interaction_matrix.csv"
    )

    plot_paths["kij_heatmap"] = _figure_path(
        output_dir, run_stamp, fig_num, "kij_interaction_heatmap"
    )
    _render_report_figure(
        plot_pairwise_matrix_heatmap,
        kij_matrix,
        labels=labels,
        title="Kij Coupling Matrix",
        colorbar_label="Normalised Kij",
        cmap=sns.color_palette("light:darkblue", as_cmap=True),
        save_path=str(plot_paths["kij_heatmap"]),
        fig_num=fig_num,
    )
    fig_numbers["kij_heatmap"] = fig_num
    fig_num += 1

    kij_terms = results.get("kij_terms", {})
    term_names = sorted(
        {
            str(term_name)
            for term_values in kij_terms.values()
            for term_name, value in term_values.items()
            if term_name != "Kij" and np.isscalar(value)
        }
    )
    for term_name in term_names:
        term_pair_values = {
            pair: float(term_values[term_name])
            for pair, term_values in kij_terms.items()
            if term_name in term_values and np.isfinite(float(term_values[term_name]))
        }
        if not term_pair_values:
            continue

        artifact_term = _safe_artifact_name(term_name)
        term_matrix = _matrix_from_pair_values(term_pair_values, n_res)
        pd.DataFrame(term_matrix, index=labels, columns=labels).to_csv(
            output_dir / f"{slug}_kij_term_{artifact_term}_matrix.csv"
        )

        plot_key = f"kij_term_{artifact_term}_heatmap"
        plot_paths[plot_key] = _figure_path(
            output_dir, run_stamp, fig_num, f"kij_term_{artifact_term}_heatmap"
        )
        _render_report_figure(
            plot_pairwise_matrix_heatmap,
            term_matrix,
            labels=labels,
            title=f"Kij Term: {term_name}",
            colorbar_label=term_name,
            cmap=sns.color_palette("light:darkblue", as_cmap=True),
            save_path=str(plot_paths[plot_key]),
            fig_num=fig_num,
        )
        fig_numbers[plot_key] = fig_num
        fig_num += 1

    dei_matrix = directional_excitation_information(results["Z_samples"])
    plot_paths["dei_heatmap"] = _figure_path(
        output_dir, run_stamp, fig_num, "directional_excitation_information"
    )
    df = pd.DataFrame(dei_matrix, columns=labels, index=labels)
    df.to_csv(output_dir / f"{slug}_directional_excitation_information.csv")
    _render_report_figure(
        plot_directional_excitation_information_heatmap,
        dei_matrix,
        labels=labels,
        save_path=str(plot_paths["dei_heatmap"]),
        fig_num=fig_num,
    )
    fig_numbers["dei_heatmap"] = fig_num
    fig_num += 1

    # Save DEI matrix as CSV for clustermap generation
    dei_csv_path = output_dir / f"{slug}_directional_excitation_information.csv"
    dei_df = pd.DataFrame(dei_matrix, index=labels, columns=labels)
    dei_df.to_csv(dei_csv_path)

    # Generate DEI clustermap from CSV
    plot_paths["dei_clustermap_png"] = _figure_path(
        output_dir, run_stamp, fig_num, "directional_excitation_information_clustermap"
    )
    plot_paths["dei_clustermap_pdf"] = plot_paths["dei_clustermap_png"].with_suffix(".pdf")
    _render_report_figure(
        plot_directional_excitation_information_clustermap,
        str(dei_csv_path),
        save_png=str(plot_paths["dei_clustermap_png"]),
        save_pdf=str(plot_paths["dei_clustermap_pdf"]),
    )
    fig_numbers["dei_clustermap_png"] = fig_num
    fig_num += 1

    mi_matrix = mutual_information_matrix(results["Z_samples"])
    mi_csv_path = output_dir / f"{slug}_mutual_information.csv"
    mi_df = pd.DataFrame(mi_matrix, index=labels, columns=labels)
    mi_df.to_csv(mi_csv_path)

    plot_paths["mi_heatmap"] = _figure_path(
        output_dir, run_stamp, fig_num, "mutual_information"
    )
    _render_report_figure(
        plot_pairwise_matrix_heatmap,
        mi_matrix,
        labels=labels,
        title="Mutual Information Matrix",
        colorbar_label="MI (bits)",
        cmap=sns.color_palette("light:darkblue", as_cmap=True),
        save_path=str(plot_paths["mi_heatmap"]),
        fig_num=fig_num,
    )
    fig_numbers["mi_heatmap"] = fig_num
    fig_num += 1

    try:
        allostery_outputs = analyse_motif_allostery(
            Z_samples=results["Z_samples"],
            residue_labels=labels,
            min_support=20,
            save_dir=str(output_dir),
            save_prefix="motif_allostery",
            pdb_path=cfg.pdb_path,
            chain=chain,
            protein_only=cfg.protein_only,
        )
        for plot_key, path_key in [
            ("dimer_affinity_heatmap", "dimer_affinity_matrix_heatmap_png"),
            ("trimer_affinity_heatmap", "trimer_affinity_matrix_heatmap_png"),
        ]:
            affinity_heatmap_path = allostery_outputs.get("saved_paths", {}).get(path_key)
            if affinity_heatmap_path:
                affinity_heatmap = Path(affinity_heatmap_path)
                if affinity_heatmap.exists():
                    plot_paths[plot_key] = affinity_heatmap
                    fig_numbers[plot_key] = fig_num
                    fig_num += 1
    except Exception as exc:
        print(f"[WARNING] Motif allostery analysis skipped: {exc}")

    joint_counts = joint_state_matrices(results["Z_samples"])
    plot_paths["joint_states"] = _figure_path(
        output_dir, run_stamp, fig_num, "joint_state_frequencies"
    )
    _render_report_figure(
        plot_joint_state_frequencies,
        matrices=list(joint_counts),
        titles=[
            "(0,0): Both ground",
            "(0,1): i ground, j excited",
            "(1,0): i excited, j ground",
            "(1,1): Both excited",
        ],
        labels=labels,
        save_path=str(plot_paths["joint_states"]),
        fig_num=fig_num,
    )
    fig_numbers["joint_states"] = fig_num
    fig_num += 1

    plot_paths["ground_comparison"] = _figure_path(output_dir, run_stamp, fig_num, "ground_comparison")
    _render_report_figure(
        plot_ground_state_comparison,
        results_df,
        results,
        molecule_name=_figure_label(fig_num),
        save_path=str(plot_paths["ground_comparison"]),
        fig_num=fig_num,
    )
    fig_numbers["ground_comparison"] = fig_num
    fig_num += 1

    return AnalysisOutputs(
        results=results,
        results_df=results_df,
        info_df=info_df,
        dei_matrix=dei_matrix,
        mi_matrix=mi_matrix,
        joint_counts=joint_counts,
        plot_paths=plot_paths,
        fig_numbers=fig_numbers,
        run_stamp=run_stamp,
    )


def _scaled_report_image(path: Path, max_width: float = 6.5 * inch) -> Image:
    reader = ImageReader(str(path))
    iw, ih = reader.getSize()
    scale = min(max_width / float(iw), 1.0)
    return Image(str(path), width=iw * scale, height=ih * scale)

@profile
def generate_report(
    output_path: Path,
    molecule_name: str,
    cfg: Config,
    chain: Optional[str],
    outputs: AnalysisOutputs,
) -> None:
    styles = getSampleStyleSheet()
    title_style = styles["Title"]
    heading = styles["Heading2"]
    body = styles["BodyText"]

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=letter,
        topMargin=0.5 * inch,
        bottomMargin=0.5 * inch,
        leftMargin=0.7 * inch,
        rightMargin=0.7 * inch,
    )
    story: list = []

    story.append(Paragraph(f"{molecule_name} qprotsim report", title_style))
    story.append(Spacer(1, 0.12 * inch))
    table_index = 1
    section_index = 1

    summary_data = [
        ["Field", "Value"],
        ["PDB path", cfg.pdb_path],
        ["Chain", chain if chain is not None else "--"],
        ["Shots", f"{cfg.shots}"],
        ["Device", cfg.pl_device],
        ["Residues", f"{len(outputs.results_df)}"],
        ["Generated", _format_run_stamp(outputs.run_stamp)],
        ["Mean energy", f"{np.mean(outputs.results['E']):.4f}"],
        ["Std energy", f"{np.std(outputs.results['E']):.4f}"],
    ]
    summary_table = Table(summary_data, colWidths=[1.6 * inch, 4.8 * inch])
    summary_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ]
        )
    )

    story.append(Paragraph(f"{section_index}. Run summary", heading))
    section_index += 1
    story.append(summary_table)
    story.append(Paragraph(f"Table {table_index}. Run summary.", body))
    table_index += 1
    story.append(Spacer(1, 0.14 * inch))

    # Add structure visualizations if available (RCSB download first, then PyMOL)
    if "structure_rcsb" in outputs.plot_paths:
        structure_path = outputs.plot_paths["structure_rcsb"]
        if structure_path and structure_path.exists():
            structure_fig_num = outputs.fig_numbers.get("structure_rcsb", 1)
            story.append(Paragraph(f"Figure {structure_fig_num}: Protein Structure (RCSB)", heading))
            story.append(_scaled_report_image(structure_path, max_width=5.0 * inch))
            story.append(
                Paragraph(
                    f"Figure {structure_fig_num}. Protein Structure (RCSB). "
                    f"Official structure image from RCSB PDB database.",
                    body,
                )
            )
            story.append(Spacer(1, 0.14 * inch))

    if "structure_pymol" in outputs.plot_paths:
        pymol_path = outputs.plot_paths["structure_pymol"]
        if pymol_path and pymol_path.exists():
            pymol_fig_num = outputs.fig_numbers.get("structure_pymol", 2)
            story.append(Paragraph(f"Figure {pymol_fig_num}: Protein Structure (PyMOL Rendering)", heading))
            story.append(_scaled_report_image(pymol_path, max_width=5.0 * inch))
            story.append(
                Paragraph(
                    f"Figure {pymol_fig_num}. Protein Structure (PyMOL Rendering). "
                    f"High-quality ray-traced rendering with spectrum coloring from N-terminus (blue) to C-terminus (red).",
                    body,
                )
            )
            story.append(Spacer(1, 0.14 * inch))

    all_info = outputs.info_df.sort_values("res_index")
    table_data = [["Residue", "KL_ground (normalised)", "KL_exc (normalised)"]]
    for _, row in all_info.iterrows():
        res_label = f"{int(row['res_index'])}-{row['resname']}"
        table_data.append([res_label, f"{row['KL_ground']:.4f}", f"{row['KL_exc']:.4f}"])
    info_table = Table(table_data, colWidths=[1.6 * inch, 2.4 * inch, 2.4 * inch])
    info_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ]
        )
    )
    story.append(Paragraph(f"{section_index}. Normalised residue KL metrics", heading))
    section_index += 1
    story.append(info_table)
    story.append(Paragraph(f"Table {table_index}. Normalised KL metrics for all residues.", body))
    table_index += 1
    story.append(Spacer(1, 0.14 * inch))

    plot_sections = [
        (
            "Energy distribution",
            "energy_distribution",
            "Histogram of total energies (kJ mol^-1) from quantum sampling; x-axis is energy, y-axis is frequency.",
        ),
        (
            "Energy density",
            "energy_density",
            "Density of sampled total energies (kJ mol^-1); x-axis is energy and y-axis is the probability density.",
        ),
        (
            "Residue KL divergence",
            "kl_divergence",
            "Scatter plot comparing normalised KL_ground against KL_exc for each residue; the dashed diagonal marks equal ground and excited contribution.",
        ),
        (
            "Residue KL divergence line comparison",
            "kl_divergence_line",
            "Per-residue normalised KL_ground and KL_exc shown as connected line traces; dotted connectors show the difference within each residue.",
        ),
        (
            "Stability vs energetic sensitivity",
            "stability_vs_sensitivity_quadrants",
            "Residue scatter plot of ground-state stability versus normalised excited-state sensitivity, split by reference lines into four quadrants.",
        ),
        (
            "Dual axis comparison",
            "dual_axis",
            "Left axis: ground-state probability; right axis: normalised excited-state sensitivity across residues.",
        ),
        (
            "Kij coupling heatmap",
            "kij_heatmap",
            "Residue-residue Kij coupling matrix from the coupling-factor model.",
        ),
        (
            "Directional excitation information heatmap",
            "dei_heatmap",
            "Residue-residue directional excitation information (bits); brighter cells indicate stronger coupling between residue pairs.",
        ),
        (
            "Directional excitation information clustermap",
            "dei_clustermap_png",
            "Hierarchical clustermap of residue-residue directional excitation information using Manhattan distance.",
        ),
        (
            "Mutual information heatmap",
            "mi_heatmap",
            "Symmetric residue-residue mutual information in bits, computed as KL divergence between joint and independent state distributions.",
        ),
        (
            "Dimer-region affinity heatmap",
            "dimer_affinity_heatmap",
            "Dimer-region symmetric affinity matrix computed from peptide-bonded residue pairs. Dimer definitions are available in motif_allostery_dimer_region_summary.csv and motif_allostery_residue_memberships.csv.",
        ),
        (
            "Trimer-region affinity heatmap",
            "trimer_affinity_heatmap",
            "Trimer-region symmetric affinity matrix computed from peptide trimers and peptide/non-covalent triangulated triads. Trimer and triad definitions are available in motif_allostery_trimer_region_summary.csv and motif_allostery_residue_memberships.csv.",
        ),
        (
            "Joint state frequencies",
            "joint_states",
            "Joint state frequencies for residue pairs; each panel is a state combination (0=ground, 1=excited).",
        ),
        (
            "Ground state comparison",
            "ground_comparison",
            "Per-residue ground-state probability before vs after entanglement shown as connected line traces on one y-axis.",
        ),
    ]

    kij_term_sections = []
    for key in outputs.plot_paths:
        if not (key.startswith("kij_term_") and key.endswith("_heatmap")):
            continue
        term_key = key[len("kij_term_") : -len("_heatmap")]
        term_label = term_key.replace("_", " ")
        kij_term_sections.append(
            (
                f"Kij term heatmap: {term_label}",
                key,
                f"Residue-residue heatmap for Kij term {term_label}.",
            )
        )

    plot_sections.extend(kij_term_sections)
    plot_sections.sort(key=lambda item: outputs.fig_numbers.get(item[1], 10**9))

    for title, key, caption in plot_sections:
        path = outputs.plot_paths.get(key)
        if not path or not path.exists():
            continue
        fig_num = outputs.fig_numbers.get(key)
        if fig_num is None:
            # Skip if no figure number assigned (shouldn't happen)
            print(f"[WARNING] No figure number for {key}, skipping...")
            continue
        story.append(
            KeepTogether(
                [
                    Paragraph(f"Figure {fig_num}: {title}", heading),
                    _scaled_report_image(path),
                    Paragraph(f"Figure {fig_num}. {title}. {caption}", body),
                ]
            )
        )
        story.append(Spacer(1, 0.14 * inch))

    story.append(Paragraph(f"{section_index}. Notes", heading))
    story.append(
        Paragraph(
            "All plots use seaborn styling for consistent academic presentation.",
            body,
        )
    )
    if "dimer_affinity_heatmap" in outputs.plot_paths or "trimer_affinity_heatmap" in outputs.plot_paths:
        story.append(
            Paragraph(
                "Dimer and trimer motif definitions are available in "
                "motif_allostery_dimer_region_summary.csv and motif_allostery_trimer_region_summary.csv; "
                "residue memberships are available in "
                "motif_allostery_residue_memberships.csv.",
                body,
            )
        )

    doc.build(story)


def run_analysis_from_yaml(
    config_path: Path,
    output_dir: Optional[Path] = None,
    report_filename: Optional[str] = None,
) -> Path:
    """
    Run analysis from YAML config file.
    
    Parameters
    ----------
    config_path : Path
        Path to YAML configuration file
    output_dir : Optional[Path]
        Override output directory (takes precedence over config file)
    report_filename : Optional[str]
        Override report filename (takes precedence over config file).
        The report will always be placed at the root of output_dir.
    
    Returns
    -------
    Path
        Path to the generated report PDF
    """
    config = load_yaml_config(config_path)
    cfg, molecule_name, chain, bins, config_output_dir, config_report_filename = _config_from_data(config)
    run_stamp = _timestamp()
    output_dir_resolved, report_path_resolved = _resolve_output_paths(
        config_path,
        output_dir,
        config_output_dir,
        report_filename,
        config_report_filename,
        molecule_name,
        run_stamp=run_stamp,
    )
    report_path_resolved.parent.mkdir(parents=True, exist_ok=True)
    
    # print(f"[INFO] Output directory: {output_dir_resolved}")
    # print(f"[INFO] Report will be saved to: {report_path_resolved}")
    # print(f"[INFO] Using chain from config: {chain if chain is not None else '--'}")

    outputs = run_analysis(
        cfg=cfg,
        molecule_name=molecule_name,
        chain=chain,
        output_dir=output_dir_resolved,
        bins=bins,
        run_stamp=run_stamp,
    )
    print("outputs computed! Generating report...")
    generate_report(
        output_path=report_path_resolved,
        molecule_name=molecule_name,
        cfg=cfg,
        chain=chain,
        outputs=outputs,
    )
    return report_path_resolved


def run_batch_analysis_from_yaml(
    config_path: Path,
    output_dir: Optional[Path] = None,
) -> dict[str, Path]:
    """
    Run batch analysis for multiple proteins from a single YAML config.
    
    A batch YAML config must contain batch_mode=true and a pdb_codes list.
    Common analysis settings apply to every downloaded structure.
    
    Example YAML format for batch processing with PDB codes:
    ```yaml
    batch_mode: true
    pdb_codes: ['1L2Y', '1AKG', '1B03']
    download_dir: /path/to/pdbs
    force_download: false
    verify_checksum: true
    
    # Common settings for all analyses
    protein_only: true
    shots: 1000000
    chain: A
    bins: 10000
    ```
    
    Parameters
    ----------
    config_path : Path
        Path to YAML configuration file
    output_dir : Optional[Path]
        Base output directory for all analyses
        
    Returns
    -------
    dict[str, Path]
        Dictionary mapping molecule names to their report paths
    """
    config = load_yaml_config(config_path)
    
    # Check if this is a batch configuration
    if not config.get("batch_mode", False):
        # Single analysis mode
        report = run_analysis_from_yaml(config_path, output_dir)
        return {"single": report}
    
    # Batch mode processing
    pdb_codes = config.get("pdb_codes", [])
    if not pdb_codes:
        raise ValueError("batch_mode is true but no pdb_codes provided")
    
    print(f"\n{'='*70}")
    print(f"BATCH PROCESSING MODE: {len(pdb_codes)} proteins")
    print(f"{'='*70}\n")
    
    # Download all PDB files first (if using PDB codes)
    from qprotsim.helper import batch_fetch_pdb
    
    download_dir = config.get("download_dir", None)
    force_download = config.get("force_download", False)
    verify_checksum = config.get("verify_checksum", True)
    
    print("[STEP 1/2] Downloading PDB files...")
    pdb_paths = batch_fetch_pdb(
        pdb_codes=pdb_codes,
        output_dir=download_dir,
        force_download=force_download,
        verify_checksum=verify_checksum,
        delay=0.5
    )
    
    print(f"\n[STEP 2/2] Running analyses for {len(pdb_paths)} proteins...")
    print("="*70)
    
    # Common settings.
    common_settings = {
        "protein_only": config.get("protein_only", True),
        "shots": config.get("shots", 100_000),
        "pl_device": config.get("pl_device", "default.qubit"),
        "chain": config.get("chain", "A"),
        "bins": config.get("bins", 1000),
    }
    
    # Determine base output directory for batch
    # Priority: function parameter > config file > default
    if output_dir is not None:
        base_output_dir = output_dir
    elif "output_dir" in config:
        base_output_dir = Path(config["output_dir"])
    else:
        base_output_dir = config_path.parent / "qprotsim_batch_outputs"
    
    # Get report filename pattern from config (supports {molecule_name} and {pdb_code})
    report_filename_pattern = config.get("report_file_name", "{pdb_code}_report.pdf")
    
    # Process each protein
    reports = {}
    for i, (pdb_code, pdb_path) in enumerate(pdb_paths.items(), 1):
        print(f"\n[{i}/{len(pdb_paths)}] Analyzing {pdb_code}...")

        # Force garbage collection
        gc.collect()
        
        try:
            run_stamp = _timestamp()
            # Create config for this protein
            molecule_name = config.get("molecule_names", {}).get(pdb_code, pdb_code)
            
            protein_config = {
                "pdb_path": str(pdb_path),
                "molecule_name": molecule_name,
                **common_settings
            }
            
            # Create a temporary config dict and process
            cfg, mol_name, chain, bins, _, _ = _config_from_data(protein_config)
            
            # Each protein gets its own subdirectory
            protein_output = base_output_dir / pdb_code.lower()
            protein_output.mkdir(parents=True, exist_ok=True)
            
            # Generate report filename from pattern
            # Support {molecule_name} and {pdb_code} placeholders
            report_filename = report_filename_pattern.format(
                molecule_name=mol_name.lower().replace(" ", "_"),
                pdb_code=pdb_code.lower(),
            )
            report_file = protein_output / report_filename
            
            # Run analysis
            outputs = run_analysis(
                cfg=cfg,
                molecule_name=mol_name,
                chain=chain,
                output_dir=protein_output,
                bins=bins,
                run_stamp=run_stamp,
            )
            
            # Generate report
            generate_report(
                output_path=report_file,
                molecule_name=mol_name,
                cfg=cfg,
                chain=chain,
                outputs=outputs,
            )
            
            reports[pdb_code] = report_file
            print(f"[SUCCESS] Generated report: {report_file}")
            
            # Force memory cleanup after each report

            # Get memory usage before cleanup
            process = psutil.Process()
            cpu_mem_before = process.memory_info().rss / (1024**3)  # GB
            
            gpu_mem_before = None
            try:
                from qprotsim.helper import HAVE_CUPY
                if HAVE_CUPY:
                    import cupy as cp
                    mempool = cp.get_default_memory_pool()
                    gpu_mem_before = mempool.used_bytes() / (1024**3)  # GB
            except:
                pass
            
            print(f"[MEMORY] Before cleanup - CPU: {cpu_mem_before:.2f} GB", end="")
            if gpu_mem_before is not None:
                print(f", GPU: {gpu_mem_before:.2f} GB")
            else:
                print()
            
            # === CRITICAL: Explicitly delete large objects to break references ===
            # 1. Delete the outputs object and its large arrays
            if 'outputs' in locals():
                # Clear large numpy arrays from results dict
                if hasattr(outputs, 'results') and outputs.results:
                    if 'Z_samples' in outputs.results:
                        del outputs.results['Z_samples']
                    if 'E' in outputs.results:
                        del outputs.results['E']
                    outputs.results.clear()
                
                # Clear dataframes
                if hasattr(outputs, 'results_df'):
                    del outputs.results_df
                if hasattr(outputs, 'info_df'):
                    del outputs.info_df
                
                # Clear matrices
                if hasattr(outputs, 'dei_matrix'):
                    del outputs.dei_matrix
                if hasattr(outputs, 'mi_matrix'):
                    del outputs.mi_matrix
                if hasattr(outputs, 'joint_counts'):
                    del outputs.joint_counts
                
                del outputs
            
            # 2. Delete config object
            if 'cfg' in locals():
                del cfg
            
            # 3. Close all matplotlib figures
            import matplotlib.pyplot as plt
            plt.close('all')
            
            # 4. Try to clear PyMOL state if it was used
            try:
                import pymol
                from pymol import cmd
                cmd.reinitialize()
            except:
                pass
            
            # 5. Force aggressive garbage collection
            # Collect each generation explicitly for thorough cleanup
            collected = gc.collect(0)  # Collect generation 0 (youngest)
            collected += gc.collect(1)  # Collect generation 1
            collected += gc.collect(2)  # Collect generation 2 (oldest, full collection)
            print(f"[MEMORY] Garbage collected {collected} objects")
            
            # 7. Clear GPU memory if available
            if gpu_mem_before is not None:
                try:
                    mempool.free_all_blocks()
                    pinned_mempool = cp.get_default_pinned_memory_pool()
                    pinned_mempool.free_all_blocks()
                except:
                    pass
            
            # Get memory usage after cleanup
            cpu_mem_after = process.memory_info().rss / (1024**3)  # GB
            gpu_mem_after = None
            if gpu_mem_before is not None:
                try:
                    gpu_mem_after = mempool.used_bytes() / (1024**3)  # GB
                except:
                    pass
            
            print(f"[MEMORY] After cleanup  - CPU: {cpu_mem_after:.2f} GB", end="")
            if gpu_mem_after is not None:
                print(f", GPU: {gpu_mem_after:.2f} GB")
                print(f"[MEMORY] Freed - CPU: {cpu_mem_before - cpu_mem_after:.2f} GB, GPU: {gpu_mem_before - gpu_mem_after:.2f} GB")
            else:
                print()
                print(f"[MEMORY] Freed - CPU: {cpu_mem_before - cpu_mem_after:.2f} GB")
            
        except Exception as e:
            print(f"[ERROR] Failed to analyze {pdb_code}: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"BATCH PROCESSING COMPLETE: {len(reports)}/{len(pdb_paths)} successful")
    print("="*70)
    
    return reports
