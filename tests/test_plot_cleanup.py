import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from qprotsim import plots
from qprotsim import report_analysis


def test_histogram_saves_and_closes_its_figure(tmp_path):
    existing_fig = plt.figure()
    try:
        output = tmp_path / "histogram.png"
        plots.plot_energy_distribution_histogram(
            np.array([-3.0, -2.0, -1.0]), "test", bins=3, save_path=str(output)
        )
        assert output.is_file()
        assert output.with_suffix(".pdf").is_file()
        assert plt.get_fignums() == [existing_fig.number]
    finally:
        plt.close(existing_fig)


def test_histogram_closes_figure_when_save_fails(monkeypatch, tmp_path):
    def fail_save(*args, **kwargs):
        raise OSError("save failed")

    monkeypatch.setattr(plots, "_save_figure_with_pdf", fail_save)
    before = plt.get_fignums()
    with pytest.raises(OSError, match="save failed"):
        plots.plot_energy_distribution_histogram(
            np.array([-3.0, -2.0, -1.0]),
            "test",
            bins=3,
            save_path=str(tmp_path / "histogram.png"),
        )
    assert plt.get_fignums() == before


@pytest.mark.parametrize(
    "plot",
    [plots.plot_energy_distribution, plots.plot_energy_density],
)
def test_energy_plots_return_figure_and_axes_after_saving(plot, monkeypatch, tmp_path):
    monkeypatch.setattr(plots, "_save_figure_with_pdf", lambda *args, **kwargs: None)
    fig, ax = plot(
        np.array([-3.0, -2.0, -1.0]),
        "test",
        bins=3,
        save_path=str(tmp_path / "energy.png"),
    )
    assert ax.figure is fig
    assert fig.number not in plt.get_fignums()


def test_scatter_annotations_do_not_close_figure():
    fig, ax = plt.subplots()
    try:
        plots._annotate_scatter_labels(ax, [0.5], [0.5], ["residue"])
        assert fig.number in plt.get_fignums()
        assert len(ax.texts) == 1
    finally:
        plt.close(fig)


def test_report_heatmap_saves_both_formats_and_keeps_existing_figures(tmp_path):
    existing_fig = plt.figure()
    try:
        baseline = tmp_path / "baseline.png"
        output = tmp_path / "heatmap.png"
        fig, ax = plots.plot_pairwise_matrix_heatmap(
            np.array([[0.0, 0.5], [0.5, 0.0]]),
            labels=["A", "B"],
            title="Heatmap",
            colorbar_label="Value",
            save_path=str(baseline),
        )
        plt.close(fig)
        del fig, ax
        report_analysis._render_report_figure(
            plots.plot_pairwise_matrix_heatmap,
            np.array([[0.0, 0.5], [0.5, 0.0]]),
            labels=["A", "B"],
            title="Heatmap",
            colorbar_label="Value",
            save_path=str(output),
        )
        assert output.is_file()
        assert output.read_bytes() == baseline.read_bytes()
        assert output.with_suffix(".pdf").is_file()
        assert plt.get_fignums() == [existing_fig.number]
    finally:
        plt.close(existing_fig)


def test_report_figure_closes_when_plot_fails():
    existing_fig = plt.figure()
    try:
        def fail_plot():
            plt.figure()
            raise OSError("plot failed")

        with pytest.raises(OSError, match="plot failed"):
            report_analysis._render_report_figure(fail_plot)
        assert plt.get_fignums() == [existing_fig.number]
    finally:
        plt.close(existing_fig)


def test_run_analysis_closes_every_report_plot(monkeypatch, tmp_path):
    cfg = SimpleNamespace(pdb_path="test.pdb", protein_only=True)
    samples = np.array([[1, -1], [-1, 1]])
    results = {
        "Z_samples": samples,
        "E": np.array([-1.0, -2.0]),
        "kijs": {(0, 1): 0.5},
        "kij_terms": {(0, 1): {"first": 0.1, "second": 0.2}},
    }
    stats = pd.DataFrame({"res_index": [1, 2], "resname": ["ALA", "GLY"]})
    info = pd.DataFrame({"res_index": [1, 2], "KL_ground": [0.5, 0.5], "KL_exc": [0.5, 0.5]})
    monkeypatch.setattr(report_analysis, "get_structure_image", lambda **kwargs: None)
    monkeypatch.setattr(report_analysis, "render_structure_pymol", lambda **kwargs: None)
    monkeypatch.setattr(report_analysis, "run_pipeline", lambda *args, **kwargs: results)
    monkeypatch.setattr(report_analysis, "export_residue_statistics", lambda *args, **kwargs: stats)
    monkeypatch.setattr(report_analysis, "export_entanglement_parameters", lambda *args, **kwargs: None)
    monkeypatch.setattr(report_analysis, "residue_information_gain", lambda *args: info)
    monkeypatch.setattr(report_analysis, "directional_excitation_information", lambda *args: np.zeros((2, 2)))
    monkeypatch.setattr(report_analysis, "mutual_information_matrix", lambda *args: np.zeros((2, 2)))
    monkeypatch.setattr(report_analysis, "joint_state_matrices", lambda *args: (np.zeros((2, 2)),) * 4)
    monkeypatch.setattr(report_analysis, "analyse_motif_allostery", lambda **kwargs: {"saved_paths": {}})

    existing_fig = plt.figure()
    plot_calls = []

    def fake_plot(*args, **kwargs):
        assert plt.get_fignums() == [existing_fig.number]
        plot_calls.append(1)
        return plt.subplots()

    for name in (
        "plot_energy_distribution_histogram",
        "plot_energy_density",
        "plot_kl_divergence_scatter",
        "plot_kl_divergence_line_comparison",
        "plot_stability_vs_sensitivity_quadrants",
        "plot_dual_axis_line_comparison",
        "plot_pairwise_matrix_heatmap",
        "plot_directional_excitation_information_heatmap",
        "plot_directional_excitation_information_clustermap",
        "plot_joint_state_frequencies",
        "plot_ground_state_comparison",
    ):
        monkeypatch.setattr(report_analysis, name, fake_plot)

    try:
        outputs = report_analysis.run_analysis(cfg, "test", None, tmp_path, bins=3)
        assert len(plot_calls) == 14
        assert {"kij_term_first_heatmap", "kij_term_second_heatmap", "dei_clustermap_png"} <= outputs.plot_paths.keys()
        assert plt.get_fignums() == [existing_fig.number]
    finally:
        plt.close(existing_fig)
