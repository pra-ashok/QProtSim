"""Verify YAML backend selection and explicit sample conventions in reports."""

from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

yaml = pytest.importorskip("yaml")
pytest.importorskip("seaborn")
pytest.importorskip("psutil")
pytest.importorskip("reportlab")
pytest.importorskip("freesasa")

from qprotsim import helper, report_analysis
from qprotsim.circuit import Config


@pytest.fixture
def local_structure(tmp_path):
    pdb_path = tmp_path / "protein.pdb"
    pdb_path.write_text("END\n", encoding="utf-8")
    return pdb_path


def write_config(tmp_path, data):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return config_path


def test_yaml_selects_qiskit_backend_and_preserves_simulator_options(tmp_path, local_structure):
    data = {
        "pdb_path": str(local_structure),
        "molecule_name": "Example protein",
        "chain": "A",
        "shots": 17,
        "backend": "qiskit",
        "qiskit_backend": "aer_simulator_statevector",
        "backend_options": {"device": "CPU", "max_parallel_threads": 2},
        "seed": 0,
    }
    parsed = report_analysis.load_yaml_config(write_config(tmp_path, data))
    cfg, name, chain, *_ = report_analysis._config_from_data(parsed)
    assert name == "Example protein"
    assert chain == "A"
    assert cfg.shots == 17
    assert cfg.backend == "qiskit"
    assert cfg.qiskit_backend == "aer_simulator_statevector"
    assert cfg.backend_options == {"device": "CPU", "max_parallel_threads": 2}
    assert cfg.seed == 0


def test_existing_yaml_keeps_pennylane_defaults(tmp_path, local_structure):
    config_path = write_config(tmp_path, {
        "pdb_path": str(local_structure),
        "molecule_name": "Example protein",
        "pl_device": "lightning.qubit",
    })
    cfg, *_ = report_analysis._config_from_data(report_analysis.load_yaml_config(config_path))
    assert cfg.backend == "pennylane"
    assert cfg.pl_device == "lightning.qubit"
    assert cfg.shots == 100_000
    assert cfg.qiskit_backend == "aer_simulator"
    assert cfg.backend_options == {}
    assert cfg.seed is None


@pytest.mark.parametrize("key,value", [
    ("backend", "invalid"),
    ("shots", 2.5),
    ("shots", "100"),
    ("shots", True),
    ("backend_options", ["method", "statevector"]),
    ("backend_options", None),
    ("seed", -1),
    ("seed", True),
])
def test_invalid_yaml_backend_values_are_rejected_without_coercion(
    tmp_path, local_structure, key, value,
):
    config_path = write_config(tmp_path, {
        "pdb_path": str(local_structure),
        "molecule_name": "Example protein",
        key: value,
    })
    with pytest.raises(ValueError):
        report_analysis._config_from_data(report_analysis.load_yaml_config(config_path))


def test_batch_yaml_forwards_backend_settings_to_every_protein(tmp_path, monkeypatch):
    pdb_paths = {}
    for code in ("1AAA", "2BBB"):
        pdb_path = tmp_path / f"{code}.pdb"
        pdb_path.write_text("END\n", encoding="utf-8")
        pdb_paths[code] = pdb_path
    batch_fetch = Mock(return_value=pdb_paths)
    monkeypatch.setattr(helper, "batch_fetch_pdb", batch_fetch)
    monkeypatch.setattr(helper, "HAVE_CUPY", False)
    monkeypatch.setattr(report_analysis.gc, "collect", lambda *args: 0)
    run_analysis = Mock(side_effect=lambda **kwargs: object())
    generate_report = Mock()
    monkeypatch.setattr(report_analysis, "run_analysis", run_analysis)
    monkeypatch.setattr(report_analysis, "generate_report", generate_report)
    config_path = write_config(tmp_path, {
        "batch_mode": True,
        "pdb_codes": list(pdb_paths),
        "backend": "qiskit",
        "qiskit_backend": "aer_simulator",
        "backend_options": {"method": "matrix_product_state", "device": "CPU"},
        "shots": 25,
        "seed": 12,
        "chain": "B",
    })

    reports = report_analysis.run_batch_analysis_from_yaml(config_path, output_dir=tmp_path / "results")

    assert set(reports) == set(pdb_paths)
    assert run_analysis.call_count == generate_report.call_count == 2
    batch_fetch.assert_called_once()
    configs = [call.kwargs["cfg"] for call in run_analysis.call_args_list]
    assert {cfg.pdb_path for cfg in configs} == {str(path) for path in pdb_paths.values()}
    for call, cfg in zip(run_analysis.call_args_list, configs):
        assert cfg.backend == "qiskit"
        assert cfg.qiskit_backend == "aer_simulator"
        assert cfg.backend_options == {"method": "matrix_product_state", "device": "CPU"}
        assert cfg.shots == 25
        assert cfg.seed == 12
        assert call.kwargs["chain"] == "B"
    assert configs[0].backend_options is not configs[1].backend_options


@pytest.mark.parametrize("z_samples", [
    np.ones((4, 2), dtype=np.int8),
    np.array([[1, -1], [-1, 1], [1, 1], [-1, -1]], dtype=np.int8),
])
def test_report_uses_binary_excitations_for_information_and_z_for_joint_states(
    tmp_path, monkeypatch, z_samples,
):
    results = {"Z_samples": z_samples.copy(), "E": np.arange(len(z_samples), dtype=float)}
    stats = pd.DataFrame({"res_index": [1, 2], "resname": ["ALA", "GLY"]})
    monkeypatch.setattr(report_analysis, "get_structure_image", lambda **kwargs: None)
    monkeypatch.setattr(report_analysis, "render_structure_pymol", lambda **kwargs: None)
    monkeypatch.setattr(report_analysis, "run_pipeline", lambda *args, **kwargs: results)
    monkeypatch.setattr(report_analysis, "export_residue_statistics", lambda *args, **kwargs: stats)
    monkeypatch.setattr(report_analysis, "export_entanglement_parameters", lambda *args, **kwargs: None)
    for name in vars(report_analysis):
        if name.startswith("plot_"):
            monkeypatch.setattr(report_analysis, name, lambda *args, **kwargs: None)

    received = {}

    def information(samples, energies):
        received["residue"] = np.array(samples, copy=True)
        return pd.DataFrame({"res_index": [1, 2], "KL_ground": [0.0, 0.0], "KL_exc": [0.0, 0.0]})

    def matrix(name):
        def calculate(samples):
            received[name] = np.array(samples, copy=True)
            return np.zeros((2, 2))
        return calculate

    def motifs(**kwargs):
        received["motif"] = np.array(kwargs["Z_samples"], copy=True)
        return {"saved_paths": {}}

    def joint(samples):
        received["joint"] = np.array(samples, copy=True)
        return tuple(np.zeros((2, 2)) for _ in range(4))

    monkeypatch.setattr(report_analysis, "residue_information_gain", information)
    monkeypatch.setattr(report_analysis, "directional_excitation_information", matrix("directional"))
    monkeypatch.setattr(report_analysis, "mutual_information_matrix", matrix("mutual"))
    monkeypatch.setattr(report_analysis, "analyse_motif_allostery", motifs)
    monkeypatch.setattr(report_analysis, "joint_state_matrices", joint)

    report_analysis.run_analysis(Config("unused.pdb"), "Example", "A", tmp_path)

    assert set(received) == {"residue", "directional", "mutual", "motif", "joint"}
    for name in ("residue", "directional", "mutual", "motif"):
        np.testing.assert_array_equal(received[name], (z_samples == -1).astype(np.uint8))
    np.testing.assert_array_equal(received["joint"], z_samples)
    np.testing.assert_array_equal(results["Z_samples"], z_samples)
    exported = pd.read_csv(tmp_path / "example_z_samples.csv", index_col="shot")
    np.testing.assert_array_equal(exported.to_numpy(), z_samples)
