import sys
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from qprotsim import plots


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
