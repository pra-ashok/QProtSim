"""
Before-and-After Equivalence Test: PennyLane 0.45.1 vs Qiskit v2.6.2
====================================================================
Run with:
    pytest -v test_qiskit_penny_equivalence.py
or directly:
    python test_qiskit_penny_equivalence.py
"""
import sys
from pathlib import Path

import numpy as np
import pennylane as qml
from qiskit.quantum_info import Statevector

# Add src directory to path to allow imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from qprotsim import circuit_pl  # Original PennyLane 0.45.1 circuit.py
from qprotsim import circuit_qk  # Ported Qiskit v2.6.2 circuit.py


def _make_test_fixture():
    """
    4-residue system designed to exercise BOTH entanglement branches:
      - Pair (0, 1): drive ~= 0.3351 (> 0.02)  -> paired CRY gates
      - Pair (1, 2): drive ~= 0.0049 (<= 0.02) -> weak IsingYY / RYY diffusion gate
      - Pair (2, 3): drive ~= 0.9960 (> 0.02)  -> paired CRY gates
      - Pair (0, 3): drive ~= 0.9980 (> 0.02)  -> paired CRY gates
    """
    resnames = ["TRP", "ILE", "PHE", "GLY"]
    E0, E1 = circuit_pl.per_residue_energies(resnames)
    # E0 = [-12.8, -10.2, -10.1, -0.1], E1 = [0.0, 0.0, 0.0, 0.0]
    exc_propensities = np.array([0.80, 0.50, 0.50, 0.10], dtype=float)
    kij_pairs = [(0, 1), (1, 2), (2, 3), (0, 3)]
    kijs = {
        (0, 1): 0.85,
        (1, 2): 0.60,
        (2, 3): 0.70,
        (0, 3): 0.40,
    }
    return exc_propensities, kij_pairs, kijs, E0, E1


def test_schedule_and_both_gate_branches_identical():
    exc, pairs, kijs, E0, E1 = _make_test_fixture()

    sched_pl = circuit_pl.compute_directed_transfer_schedule(exc, pairs, kijs, E0, E1)
    sched_qk = circuit_qk.compute_directed_transfer_schedule(exc, pairs, kijs, E0, E1)

    assert len(sched_pl) == len(sched_qk) > 0
    for a, b in zip(sched_pl, sched_qk):
        assert a.keys() == b.keys()
        for k in a:
            assert np.isclose(a[k], b[k], atol=1e-15), f"Mismatch in schedule key {k}"

    # Guarantee both CRY and YY branches are covered by this test fixture
    drives = [entry["drive"] for entry in sched_pl]
    assert any(0.0 < d <= 0.02 for d in drives), "Fixture must trigger YY diffusion branch"
    assert any(d > 0.02 for d in drives), "Fixture must trigger paired CRY branch"


def test_exact_statevector_and_probabilities_match():
    exc, pairs, kijs, E0, E1 = _make_test_fixture()
    n_qubits = len(exc)

    # 1) PennyLane exact statevector (big-endian: |q0 q1 ... q_{n-1}>)
    dev = qml.device("default.qubit", wires=n_qubits)

    @qml.qnode(dev)
    def pl_state_qnode():
        circuit_pl.prepare_state_propensity(exc, pairs, kijs, E0, E1)
        return qml.state()

    psi_pl = np.asarray(pl_state_qnode(), dtype=np.complex128)
    probs_pl = np.abs(psi_pl) ** 2

    # 2) Qiskit exact statevector (reverse_qargs aligns little-endian -> big-endian)
    qc, _ = circuit_qk.prepare_state_propensity(exc, pairs, kijs, E0, E1, return_schedule=True)
    sv_qk_big_endian = Statevector.from_instruction(qc).reverse_qargs()
    psi_qk = np.asarray(sv_qk_big_endian.data, dtype=np.complex128)
    probs_qk = sv_qk_big_endian.probabilities()

    # Exact amplitude equality (zero global phase discrepancy)
    np.testing.assert_allclose(psi_qk, psi_pl, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(probs_qk, probs_pl, atol=1e-12, rtol=0.0)


def test_deterministic_endianness_and_pauli_z_sign():
    """
    Deterministic basis state |1, 0, 1, 0> (qubits 0 and 2 excited, 1 and 3 ground).
    Verifies column i == qubit i and |0> -> +1, |1> -> -1 with zero shot noise.
    """
    exc = np.array([1.0, 0.0, 1.0, 0.0])
    E0 = np.array([-10.0, -5.0, -3.0, -1.0])
    E1 = np.zeros(4)

    cfg_pl = circuit_pl.Config(pdb_path="dummy.pdb", shots=128, pl_device="default.qubit")
    cfg_qk = circuit_qk.Config(pdb_path="dummy.pdb", shots=128, backend="statevector_sampler", seed=42)

    z_pl = circuit_pl.sample_bits(exc, [], {}, E0, E1, cfg_pl)
    z_qk = circuit_qk.sample_bits(exc, [], {}, E0, E1, cfg_qk)

    expected_row = np.array([-1, +1, -1, +1])
    assert z_pl.shape == z_qk.shape == (128, 4)
    assert np.all(z_pl == expected_row)
    assert np.all(z_qk == expected_row)

    e_pl = circuit_pl.energy_from_samples(z_pl, E0, E1)
    e_qk = circuit_qk.energy_from_samples(z_qk, E0, E1)
    # Qubits 1 (-5.0) and 3 (-1.0) are in |0> (+1); qubits 0 and 2 are in |1> (0.0)
    np.testing.assert_allclose(e_pl, -6.0, atol=1e-12)
    np.testing.assert_allclose(e_qk, -6.0, atol=1e-12)


def test_end_to_end_sampling_and_energy_distribution():
    exc, pairs, kijs, E0, E1 = _make_test_fixture()
    n_qubits = len(exc)
    shots = 100_000

    # Analytic reference from exact Qiskit statevector (already proven identical to PennyLane)
    qc, _ = circuit_qk.prepare_state_propensity(exc, pairs, kijs, E0, E1, return_schedule=True)
    probs_exact = Statevector.from_instruction(qc).reverse_qargs().probabilities()

    # Exact bitstring table (big-endian rows matching probs_exact)
    bit_table = ((np.arange(2**n_qubits)[:, None] >> np.arange(n_qubits - 1, -1, -1)) & 1)
    z_table = np.where(bit_table == 1, -1.0, 1.0)
    exact_z_exp = probs_exact @ z_table
    exact_energies = circuit_qk.energy_from_samples(z_table, E0, E1)
    exact_mean_E = float(probs_exact @ exact_energies)
    exact_std_E = float(np.sqrt(probs_exact @ ((exact_energies - exact_mean_E) ** 2)))

    # Run PennyLane (before) and Qiskit (after)
    cfg_pl = circuit_pl.Config(pdb_path="dummy.pdb", shots=shots, pl_device="default.qubit")
    cfg_qk = circuit_qk.Config(pdb_path="dummy.pdb", shots=shots, pl_device="default.qubit", seed=12345)

    z_pl, sched_pl = circuit_pl.sample_bits(exc, pairs, kijs, E0, E1, cfg_pl, return_schedule=True)
    z_qk, sched_qk = circuit_qk.sample_bits(exc, pairs, kijs, E0, E1, cfg_qk, return_schedule=True)

    assert z_pl.shape == z_qk.shape == (shots, n_qubits)
    assert set(np.unique(z_qk)).issubset({-1, 1})
    assert len(sched_pl) == len(sched_qk)

    # Per-qubit <Z_i> tolerance: 1/sqrt(100_000) ~= 0.0032 -> 0.02 is >6 sigma
    np.testing.assert_allclose(z_pl.mean(axis=0), exact_z_exp, atol=0.02)
    np.testing.assert_allclose(z_qk.mean(axis=0), exact_z_exp, atol=0.02)
    np.testing.assert_allclose(z_qk.mean(axis=0), z_pl.mean(axis=0), atol=0.025)

    # Total configuration energies
    e_pl = circuit_pl.energy_from_samples(z_pl, E0, E1)
    e_qk = circuit_qk.energy_from_samples(z_qk, E0, E1)

    # Standard error of mean energy is exact_std_E / sqrt(shots)
    se_mean = exact_std_E / np.sqrt(shots)
    assert abs(e_pl.mean() - exact_mean_E) < 6.0 * se_mean
    assert abs(e_qk.mean() - exact_mean_E) < 6.0 * se_mean
    assert abs(e_qk.mean() - e_pl.mean()) < 8.0 * se_mean
    assert abs(e_qk.std() - exact_std_E) < 0.15


def test_legacy_config_and_gpu_fallback():
    exc, pairs, kijs, E0, E1 = _make_test_fixture()

    # Positional Config invocation identical to original PennyLane callers
    cfg_legacy = circuit_qk.Config("dummy.pdb", True, 512, "lightning.gpu")
    samples = circuit_qk.sample_bits(exc, pairs, kijs, E0, E1, cfg_legacy)
    assert samples.shape == (512, len(exc))
    assert set(np.unique(samples)).issubset({-1, 1})


if __name__ == "__main__":
    test_schedule_and_both_gate_branches_identical()
    test_exact_statevector_and_probabilities_match()
    test_deterministic_endianness_and_pauli_z_sign()
    test_end_to_end_sampling_and_energy_distribution()
    test_legacy_config_and_gpu_fallback()
    print("ALL PENNYLANE -> QISKIT v2.6.2 EQUIVALENCE TESTS PASSED.")