import csv
from pathlib import Path
from typing import Optional

import numpy as np
import MDAnalysis as mda
import pandas as pd

from qprotsim.circuit_qk import (
    Config,
    energy_from_samples,
    excitation_propensities,
    per_residue_energies,
    sample_bits,
)
from qprotsim.coupling_factor import build_Kij_pairs
from qprotsim.helper import map_sasa_to_mda_resorder, per_residue_sasa
from qprotsim.plots import plot_energy_density, plot_energy_distribution_histogram


def _default_density_plot_path(molecule_name: str, plot_path: Optional[str]) -> str:
    if plot_path:
        plot_file = Path(plot_path)
        return str(plot_file.with_name(f"{plot_file.stem}_density.png"))

    molecule_slug = molecule_name.strip().lower().replace(" ", "_")
    return f"{molecule_slug}_energy_density.png"


def run_pipeline(
    cfg: Config,
    molecule_name: str,
    chain: Optional[str] = "A",
    plot_energy: bool = True,
    bins: int = 1,
    plot_path: Optional[str] = None,
    save_scatter_plot: bool = True,
    scatter_plot_path: Optional[str] = None,
):
    u = mda.Universe(cfg.pdb_path)
    sasa_dict = per_residue_sasa(cfg.pdb_path, chain=chain, protein_only=cfg.protein_only)
    alpha_raw, resnames = map_sasa_to_mda_resorder(
        u,
        sasa_dict,
        protein_only=cfg.protein_only,
        chain=chain,
    )

    print("SASA distionary alpha raw", alpha_raw)

    # excitation propensity (with auto-calibrated b0 if requested)
    p_exc = excitation_propensities(sasa_dict)
    p_exc_array, _ = map_sasa_to_mda_resorder(
        u,
        p_exc,
        protein_only=cfg.protein_only,
        chain=chain,
    )
    n_res = len(resnames)
    print("Number of residues:", n_res, len(alpha_raw))
    print("resnames:", resnames)

    df = pd.DataFrame()
    df["p_exc"] = p_exc_array

    kij_pairs, kijs, kij_terms  = build_Kij_pairs(
        u,
        protein_only=cfg.protein_only,
        chain=chain,
    )

    # energies
    E0, E1 = per_residue_energies(resnames)
    print("E0",E0)
    print("E1",E1)

    Z_samples, entanglement_schedule = sample_bits(
        df["p_exc"],
        kij_pairs,
        kijs,
        E0,
        E1,
        cfg,
        return_schedule=True,
    )
    E_total = energy_from_samples(Z_samples, E0, E1)

    if plot_energy:
        plot_energy_distribution_histogram(
            E_total,
            molecule_name=molecule_name,
            bins=bins,
            save_path=plot_path,
        )
        if save_scatter_plot:
            scatter_save_path = scatter_plot_path or _default_density_plot_path(
                molecule_name,
                plot_path,
            )
            plot_energy_density(
                E_total,
                molecule_name=molecule_name,
                bins=bins,
                save_path=scatter_save_path,
            )

    p_ground = np.mean(Z_samples == 1, axis=0)


    return {
        "E": E_total,
        "E0": E0, "E1": E1,
        "p_exc": df["p_exc"], "p_ground": p_ground, "pairs": kij_pairs, "kijs":kijs, "kij_terms": kij_terms ,
        "resnames": resnames, "alpha_SASA": alpha_raw,
        "Z_samples": Z_samples,
        "entanglement_schedule": entanglement_schedule,
    }



def export_residue_statistics(results: dict, out_csv: str = "residue_propensity.csv"):
    """Save per-residue excitation propensities and features to CSV."""
    resnames = results["resnames"]
    alpha = results["alpha_SASA"]
    p_exc = results["p_exc"]
    p_ground = results["p_ground"]

    rows = []
    for i, rname in enumerate(resnames):
        rows.append({
            "res_index": i + 1,  # 1-based indexing for PyMOL residue numbering
            "resname": rname,
            "SASA_norm": float(alpha[i]),
            "p_exc": float(p_exc[i]),
            "p_ground": float(p_ground[i])
        })

    results_df = pd.DataFrame(rows)
    results_df.sort_values(by=["res_index"], inplace=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"[saved] residue-level stats → {out_csv}")

    return results_df


def export_entanglement_parameters(
    results: dict,
    out_csv: str = "entanglement_parameters.csv",
    yy_cap: float | None = np.pi / 6.0,
    yy_probability_threshold: float = 0.02,
    apply_yy_diffusion: bool = True,
):
    """Save circuit entanglement schedule parameters to CSV."""
    schedule = results.get("entanglement_schedule", [])
    resnames = list(results.get("resnames", []))
    p_exc = np.asarray(results.get("p_exc", []), dtype=float)
    E0 = np.asarray(results.get("E0", []), dtype=float)
    E1 = np.asarray(results.get("E1", []), dtype=float)

    rows = []
    for entry in schedule:
        i = int(entry["i"])
        j = int(entry["j"])
        drive = float(entry.get("drive", 0.0))
        use_yy = bool(apply_yy_diffusion and drive <= yy_probability_threshold)
        yy_theta = np.nan
        if use_yy:
            yy_theta = 2.0 * np.arcsin(np.sqrt(np.clip(drive, 0.0, 1.0)))
            if yy_cap is not None:
                yy_theta = float(np.clip(yy_theta, -abs(yy_cap), abs(yy_cap)))

        source_resname = resnames[i] if i < len(resnames) else ""
        target_resname = resnames[j] if j < len(resnames) else ""
        rows.append({
            "layer": int(entry.get("layer", len(rows) + 1)),
            "gate_type": "IsingYY" if use_yy else "CRY_pair",
            "source_qubit": i,
            "target_qubit": j,
            "source_res_index": i + 1,
            "target_res_index": j + 1,
            "source_resname": source_resname,
            "target_resname": target_resname,
            "source_label": f"{i + 1}-{source_resname}" if source_resname else str(i + 1),
            "target_label": f"{j + 1}-{target_resname}" if target_resname else str(j + 1),
            "source_p_exc": float(p_exc[i]) if i < len(p_exc) else np.nan,
            "target_p_exc": float(p_exc[j]) if j < len(p_exc) else np.nan,
            "source_E0": float(E0[i]) if i < len(E0) else np.nan,
            "source_E1": float(E1[i]) if i < len(E1) else np.nan,
            "target_E0": float(E0[j]) if j < len(E0) else np.nan,
            "target_E1": float(E1[j]) if j < len(E1) else np.nan,
            "Kij": float(entry.get("Kij", np.nan)),
            "mu_i": float(entry.get("mu_i", np.nan)),
            "mu_j": float(entry.get("mu_j", np.nan)),
            "drive": drive,
            "denom_i": float(entry.get("denom_i", np.nan)),
            "deltaP": float(entry.get("deltaP", np.nan)),
            "discount_power": int(entry.get("discount_power", 0)),
            "discount": float(entry.get("discount", np.nan)),
            "cry_phi": float(entry.get("phi", np.nan)),
            "yy_theta": float(yy_theta) if np.isfinite(yy_theta) else np.nan,
            "apply_yy_diffusion": bool(apply_yy_diffusion),
            "yy_probability_threshold": float(yy_probability_threshold),
            "yy_cap": float(yy_cap) if yy_cap is not None else np.nan,
        })

    columns = [
        "layer",
        "gate_type",
        "source_qubit",
        "target_qubit",
        "source_res_index",
        "target_res_index",
        "source_resname",
        "target_resname",
        "source_label",
        "target_label",
        "source_p_exc",
        "target_p_exc",
        "source_E0",
        "source_E1",
        "target_E0",
        "target_E1",
        "Kij",
        "mu_i",
        "mu_j",
        "drive",
        "denom_i",
        "deltaP",
        "discount_power",
        "discount",
        "cry_phi",
        "yy_theta",
        "apply_yy_diffusion",
        "yy_probability_threshold",
        "yy_cap",
    ]
    schedule_df = pd.DataFrame(rows, columns=columns)
    schedule_df.to_csv(out_csv, index=False)
    print(f"[saved] circuit entanglement parameters -> {out_csv}")
    return schedule_df
