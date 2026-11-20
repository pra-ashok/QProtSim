"""Check native circuit equivalence and the public sample encoding contract."""

import builtins

import numpy as np
import pennylane as qml
import pytest

from qprotsim.circuit import (
    Config,
    bits_to_z_samples,
    build_qiskit_circuit,
    compute_directed_transfer_schedule,
    energy_from_samples,
    prepare_state_propensity,
    qiskit_memory_to_z_samples,
    sample_bits,
)


@pytest.fixture
def qiskit():
    return pytest.importorskip("qiskit")


@pytest.fixture
def aer(qiskit):
    return pytest.importorskip("qiskit_aer")


def test_bits_and_kets_have_the_expected_z_eigenvalues_and_energies():
    # Columns follow residue order. |0> has Z=+1 and |1> has Z=-1.
    bits = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
    samples = bits_to_z_samples(bits)
    np.testing.assert_array_equal(samples, [[1, 1], [1, -1], [-1, 1], [-1, -1]])
    np.testing.assert_array_equal(
        energy_from_samples(samples, np.array([-10.0, -4.0]), np.array([1.0, 2.0])),
        [-14.0, -8.0, -3.0, 3.0],
    )
    np.testing.assert_array_equal(bits, [[0, 0], [0, 1], [1, 0], [1, 1]])


@pytest.mark.parametrize("bits", [[0, 1], [[[0, 1]]], [[-1, 1]], [[0.5, 1]], [[np.nan, 0]]])
def test_binary_conversion_rejects_wrong_shape_or_nonbinary_values(bits):
    with pytest.raises(ValueError):
        bits_to_z_samples(bits)


def test_qiskit_memory_order_is_reversed_into_residue_order():
    # Qiskit prints q2 q1 q0; QProtSim columns must be q0 q1 q2.
    np.testing.assert_array_equal(
        qiskit_memory_to_z_samples(["001", "110", "1 00"], 3),
        [[-1, 1, 1], [1, -1, -1], [1, 1, -1]],
    )


@pytest.mark.parametrize("memory", [["01"], ["0001"], ["01x"], ["0x1"]])
def test_qiskit_memory_rejects_invalid_width_or_symbols(memory):
    with pytest.raises(ValueError):
        qiskit_memory_to_z_samples(memory, 3)


def test_energy_requires_z_eigenvalues_rather_than_guessing_binary_encoding():
    with pytest.raises(ValueError):
        energy_from_samples(np.array([[0, 1]]), np.array([-4.0, -2.0]), np.zeros(2))
    np.testing.assert_array_equal(
        energy_from_samples(np.ones((1, 2)), np.array([-4.0, -2.0]), np.zeros(2)),
        [-6.0],
    )


@pytest.mark.parametrize("samples,e0,e1", [
    ([1, -1], [-2, -3], [0, 0]),
    ([[1, -1]], [-2], [0, 0]),
    ([[1, -1]], [-2, -3], [0]),
    ([[1, np.nan]], [-2, -3], [0, 0]),
])
def test_energy_rejects_invalid_sample_and_energy_dimensions(samples, e0, e1):
    with pytest.raises(ValueError):
        energy_from_samples(np.asarray(samples), np.asarray(e0), np.asarray(e1))


@pytest.mark.parametrize("kwargs", [
    {"backend": "unknown"},
    {"backend": None},
    {"shots": 0},
    {"shots": -1},
    {"shots": 1.5},
    {"shots": True},
    {"seed": -1},
    {"seed": 1.5},
    {"seed": True},
    {"backend_options": []},
    {"backend_options": None},
    {"qiskit_backend": ""},
])
def test_invalid_backend_configuration_is_rejected(kwargs):
    with pytest.raises(ValueError):
        Config(pdb_path="unused.pdb", **kwargs)


def test_backend_option_defaults_are_independent():
    first, second = Config("first.pdb"), Config("second.pdb")
    first.backend_options["method"] = "statevector"
    assert second.backend_options == {}
    assert second.backend == "pennylane"
    assert second.pl_device == "default.qubit"


@pytest.mark.parametrize("apply_yy_diffusion,yy_cap", [(True, np.pi / 6), (True, 0.07), (False, None)])
def test_qiskit_matches_pennylane_statevector_with_ordered_transfer_gates(
    qiskit, apply_yy_diffusion, yy_cap,
):
    from qiskit.quantum_info import Statevector

    p = np.array([0.17, 0.41, 0.76])
    schedule = [
        {"i": 2, "j": 0, "phi": 0.57, "drive": 0.4},
        {"i": 1, "j": 2, "phi": 0.29, "drive": 0.015},
        {"i": 0, "j": 1, "phi": 0.81, "drive": 0.3},
    ]
    kwargs = dict(schedule=schedule, apply_yy_diffusion=apply_yy_diffusion, yy_cap=yy_cap)
    circuit = build_qiskit_circuit(p, [], {}, -np.ones(3), np.zeros(3), measure=False, **kwargs)

    @qml.qnode(qml.device("default.qubit", wires=3))
    def reference():
        prepare_state_propensity(p, [], {}, -np.ones(3), np.zeros(3), **kwargs)
        return qml.state()

    # PL amplitude axes are q0,q1,q2; Qiskit amplitude axes are q2,q1,q0.
    actual = Statevector.from_instruction(circuit).data.reshape(2, 2, 2).transpose(2, 1, 0).reshape(-1)
    np.testing.assert_allclose(actual, reference(), atol=1e-12)


def test_automatic_transfer_schedule_matches_explicit_qiskit_schedule(qiskit):
    from qiskit.quantum_info import Statevector

    p = np.array([0.18, 0.42, 0.73])
    pairs = [(0, 2), (0, 1), (1, 2)]
    kijs = {(0, 2): 0.5, (0, 1): 0.2, (1, 2): 0.3}
    e0, e1 = np.array([-2.0, -5.0, -3.0]), np.zeros(3)
    schedule = compute_directed_transfer_schedule(p, pairs, kijs, e0, e1)
    assert schedule
    automatic = build_qiskit_circuit(p, pairs, kijs, e0, e1, measure=False)
    explicit = build_qiskit_circuit(p, pairs, kijs, e0, e1, measure=False, schedule=schedule)
    np.testing.assert_allclose(Statevector.from_instruction(automatic).data,
                               Statevector.from_instruction(explicit).data, atol=1e-12)


@pytest.mark.parametrize("backend", ["pennylane", "qiskit"])
@pytest.mark.parametrize("probabilities,shots", [([1, 0, 0], 7), ([0, 0, 1], 1), ([0], 1), ([1], 5)])
def test_sampling_shape_and_deterministic_kets_match_both_backends(backend, probabilities, shots, request):
    if backend == "qiskit":
        request.getfixturevalue("aer")
    p = np.asarray(probabilities)
    cfg = Config("unused.pdb", backend=backend, shots=shots, seed=19)
    samples, schedule = sample_bits(p, [], {}, -np.ones(len(p)), np.zeros(len(p)), cfg, return_schedule=True)
    assert samples.shape == (shots, len(p))
    np.testing.assert_array_equal(samples, np.tile(1 - 2 * p, (shots, 1)))
    assert schedule == []


def test_qiskit_sampling_respects_probabilities_backend_options_and_seed(aer):
    p = np.array([0.2, 0.75])
    cfg = Config("unused.pdb", backend="qiskit", shots=6000, seed=123,
                 backend_options={"method": "matrix_product_state", "max_parallel_threads": 1})
    first = sample_bits(p, [], {}, -np.ones(2), np.zeros(2), cfg)
    second = sample_bits(p, [], {}, -np.ones(2), np.zeros(2), cfg)
    np.testing.assert_array_equal(first, second)
    np.testing.assert_allclose((first == -1).mean(axis=0), p, atol=0.025)
    assert set(np.unique(first)) == {-1, 1}


def test_qiskit_accepts_numpy_integer_shots_and_seed(aer):
    cfg = Config("unused.pdb", backend="qiskit", shots=np.int64(2), seed=np.int64(42))
    samples = sample_bits(np.array([1.0]), [], {}, np.array([-2.0]), np.zeros(1), cfg)
    np.testing.assert_array_equal(samples, [[-1], [-1]])


def test_qiskit_backend_errors_do_not_silently_run_pennylane(aer):
    cfg = Config("unused.pdb", backend="qiskit", qiskit_backend="does_not_exist", shots=2)
    with pytest.raises(Exception, match="does_not_exist|backend"):
        sample_bits(np.array([0.5]), [], {}, np.array([-2.0]), np.zeros(1), cfg)


def test_invalid_qiskit_options_are_reported_instead_of_ignored(aer):
    cfg = Config("unused.pdb", backend="qiskit", shots=2,
                 backend_options={"not_an_aer_option": True})
    with pytest.raises(ValueError, match="not_an_aer_option"):
        sample_bits(np.array([0.5]), [], {}, np.array([-2.0]), np.zeros(1), cfg)


def test_pennylane_can_sample_when_qiskit_imports_are_unavailable(monkeypatch):
    original_import = builtins.__import__

    def without_qiskit(name, *args, **kwargs):
        if name.split(".")[0] in {"qiskit", "qiskit_aer"}:
            raise ModuleNotFoundError(f"No module named '{name}'", name=name)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_qiskit)
    cfg = Config("unused.pdb", shots=2)
    np.testing.assert_array_equal(
        sample_bits(np.array([0.0]), [], {}, np.array([-2.0]), np.zeros(1), cfg),
        [[1], [1]],
    )


def test_qiskit_missing_dependency_has_installation_hint(monkeypatch):
    original_import = builtins.__import__

    def without_qiskit(name, *args, **kwargs):
        if name == "qiskit" or name.startswith("qiskit."):
            raise ModuleNotFoundError("No module named 'qiskit'", name="qiskit")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_qiskit)
    with pytest.raises(ImportError, match=r"\[qiskit\]"):
        build_qiskit_circuit(np.array([0.5]), [], {}, np.array([-2.0]), np.zeros(1))
