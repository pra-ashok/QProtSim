# Porting `circuit.py` from PennyLane 0.45.1 to Qiskit v2.6.2

`circuit.py` can be ported cleanly from **PennyLane 0.45.1** to **Qiskit v2.6.2** because all classical MDAnalysis, hydrophobicity, Boltzmann propensity, schedule computation, and energy accumulation logic is framework-agnostic. Only three quantum-specific functions (`make_device`, `prepare_state_propensity`, and `sample_bits`) and the `Config` execution settings require conversion.

---

## 1. Conversion Logic & Mathematical Equivalence

### PennyLane 0.45.1 $\rightarrow$ Qiskit v2.6.2 Mapping

| Component in `circuit.py` | PennyLane 0.45.1 | Qiskit v2.6.2 Equivalent | Equivalence / Critical Note |
| :--- | :--- | :--- | :--- |
| **Single-qubit init** | `qml.RY(theta, wires=i)` | `qc.ry(theta, i)` | Exact match: $R_y(\theta) = \exp(-i\frac{\theta}{2}Y)$. |
| **Directed transfer** | `qml.CRY(phi, wires=[i, j])` | `qc.cry(phi, i, j)` | Exact physical action: `i` is control, `j` is target. |
| **Weak YY diffusion** | `qml.IsingYY(theta_yy, wires=[i, j])` | `qc.ryy(theta_yy, i, j)` | Exact match: $\exp(-i\frac{\theta}{2} Y \otimes Y)$ (same sign, factor of $\frac{1}{2}$, and qubit symmetry). |
| **Circuit container** | Implicit tape inside `@qml.qnode` | Explicit `QuantumCircuit(n_qubits)` | Qiskit has no global tape; `prepare_state_propensity` returns `qc` (or `(qc, schedule)`). |
| **Default CPU simulator** | `qml.device("default.qubit", wires=n)` + `@qml.set_shots(cfg.shots)` | `StatevectorSampler(default_shots=cfg.shots, seed=cfg.seed)` | Pure Qiskit 2.6.2 SamplerV2 primitive; no extra plugin required. |
| **Optional GPU / Aer fallback** | `lightning.gpu` $\rightarrow$ fallback to `default.qubit` | `BackendSamplerV2(backend=AerSimulator(method="statevector", device="GPU"))` $\rightarrow$ fallback to `StatevectorSampler` | Preserves `pl_device="lightning.gpu"` compatibility while adding native Qiskit backend names. |
| **Pauli-Z sampling** | `[qml.sample(qml.PauliZ(i)) for i in range(n)]` stacked with `axis=1` | `qc.measure(range(n), meas)` $\rightarrow$ `bits = pub_res.data.meas.to_bool_array(order="little")` $\rightarrow$ `np.where(bits, -1, 1)` | Preserves shape `(shots, n_qubits)`, column `i` = qubit `i`, and eigenvalue mapping $|0\rangle \rightarrow +1$, $|1\rangle \rightarrow -1$. |

### Key Technical Details (Preventing Silent Bugs)

1. **Gate Matrix & Wire-Order Verification**
   - **`RY` (`qml.RY` vs `QuantumCircuit.ry`)**: Both define
     $$R_y(\theta) = \begin{pmatrix} \cos(\theta/2) & -\sin(\theta/2) \\ \sin(\theta/2) & \cos(\theta/2) \end{pmatrix}$$
     For $\theta_i = 2\arcsin(\sqrt{p_i})$, both map $|0\rangle \mapsto \sqrt{1-p_i}|0\rangle + \sqrt{p_i}|1\rangle$.
   - **`CRY` (`qml.CRY` vs `QuantumCircuit.cry`)**: In PennyLane, `wires=[i, j]` sets wire `i` as control and wire `j` as target (`CRY` applies $R_y(\phi)$ to `j` conditioned on `i` being $|1\rangle$). In Qiskit, `qc.cry(phi, i, j)` takes `(theta, control_qubit, target_qubit)`. Although Qiskit's textbook $4 \times 4$ matrix uses little-endian tensor indexing ($q_1 \otimes q_0$) while PennyLane uses big-endian ($q_0 \otimes q_1$), the **physical unitary action on logical qubit indices `(i, j)` is strictly identical** with zero global phase difference.
   - **`IsingYY` vs `RYY` (`qml.IsingYY` vs `QuantumCircuit.ryy`)**: Both implement $\exp(-i\frac{\theta}{2} Y \otimes Y)$ with diagonal $\cos(\theta/2)$ and anti-diagonal $(+i\sin(\theta/2), -i\sin(\theta/2), -i\sin(\theta/2), +i\sin(\theta/2))$. No angle negation or factor-of-2 adjustment is needed.

2. **Endianness & Pauli-Z Eigenvalue Conversion (Most Common Porting Pitfall)**
   - In `circuit.py`, `sample_bits` returns an array `Z_samples` of shape `(shots, n_qubits)` where column `i` corresponds to residue/qubit `i`, and `energy_from_samples` computes `ground = (Z_samples == +1)`.
   - In Qiskit v2.6.2, `BitArray.to_bool_array(order="big")` (the **default**) places the most-significant bit (qubit `n - 1`) at column `0`, which would silently reverse the residue ordering.
   - Passing **`order="little"`** (`pub_result.data.meas.to_bool_array(order="little")`) places classical bit `0` (qubit `0`) at column `0` and classical bit `i` (qubit `i`) at column `i`.
   - Converting boolean measurement outcomes (`False` for $|0\rangle$, `True` for $|1\rangle$) via `np.where(bits, -1, 1)` reproduces PennyLane's exact `qml.sample(qml.PauliZ(i))` convention ($|0\rangle \rightarrow +1$, $|1\rangle \rightarrow -1$) and guarantees `energy_from_samples` produces identical energies.

3. **Zero-Breakage API Compatibility**
   - `Config` preserves all four original fields in their exact positional order (`pdb_path`, `protein_only`, `shots`, `pl_device`) and adds optional `backend: str | None = None` and `seed: int | None = None`. Any existing caller passing `Config(pdb_path=..., shots=..., pl_device="default.qubit")` or `"lightning.gpu"` continues working without modification.
   - **One unavoidable signature change**: Because Qiskit has no implicit tape queuing, `prepare_state_propensity` accepts an optional `qc: QuantumCircuit | None = None` and returns `(qc, schedule)` when `return_schedule=True` or `qc` when `return_schedule=False` (in PennyLane it returned `schedule` or `None` while mutating the active QNode tape). `sample_bits` keeps its **exact** original signature and return contract.

---

## 2. Official API References

### Qiskit v2.6.2 SDK
- `StatevectorSampler`: <https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.primitives.StatevectorSampler>
- `BackendSamplerV2`: <https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.primitives.BackendSamplerV2>
- `BitArray.to_bool_array`: <https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.primitives.BitArray>
- `RYGate` / `CRYGate` / `RYYGate`: <https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.circuit.library.RYYGate>
- `Statevector.reverse_qargs`: <https://quantum.cloud.ibm.com/docs/api/qiskit/qiskit.quantum_info.Statevector>

### PennyLane 0.45.1
- `qml.CRY`: <https://docs.pennylane.ai/en/stable/code/api/pennylane.CRY.html>
- `qml.IsingYY`: <https://docs.pennylane.ai/en/stable/code/api/pennylane.IsingYY.html>

---

## 3. Ported `circuit.py` (Qiskit v2.6.2)

All non-quantum utilities from the original `circuit.py` are preserved verbatim so downstream imports never break.

```python
"""
Quantum Energy Propensity for Proteins -- KE-aware & ground-skewed (Qiskit v2.6.2)
==================================================================================

- Qubit residue r: |0> = stable solvation, |1> = excited solvation.
- p_exc(r) = normalised(SASA_i) or 1-normalised(SASA_i) -> depending on hydrophobicity
- E0(r) negative; E1(r) = E0(r) + gap(r). The gap depends on the favourable solvation energy.
- Residue-residue Kij coupling factors drive directed excitation-transfer layers.
- Strong directed transfers use paired CRY gates; weak transfers can use YY diffusion.

Author: Pratik (Ported to Qiskit v2.6.2)
"""
from dataclasses import dataclass
import sys
from typing import Dict, List, Sequence, Tuple, Mapping, Any

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import distance_array
from qiskit import ClassicalRegister, QuantumCircuit
from qiskit.primitives import BackendSamplerV2, BaseSamplerV2, StatevectorSampler
from qiskit.transpiler import generate_preset_pass_manager

from qprotsim.info_theory import debug_mode


@dataclass
class Config:
    """Structure and sampling settings, using CPU execution by default."""

    pdb_path: str
    protein_only: bool = True

    # shots and device/backend
    shots: int = 100_000
    # Kept in position 4 for backward compatibility with legacy callers.
    # Accepts legacy names ("default.qubit", "lightning.gpu") as well as
    # Qiskit names ("statevector_sampler", "aer_cpu", "aer_gpu").
    pl_device: str = "statevector_sampler"
    backend: str | None = None
    seed: int | None = None


# --------------------- Residue properties -------------------
# Fauchere and Pliska, 1983
# Energies wrt bulk octanol
energies_bulk_octabol = {
    "F": {"e0": -10.1, "e1": 0.0},
    "I": {"e0": -10.2, "e1": 0.0},
    "W": {"e0": -12.8, "e1": 0.0},
    "L": {"e0": -9.64, "e1": 0.0},
    "Y": {"e0": -5.44, "e1": 0.0},
    "M": {"e0": -6.97, "e1": 0.0},
    "V": {"e0": -6.91, "e1": 0.0},
    "P": {"e0": -4.08, "e1": 0.0},
    "C": {"e0": -8.73, "e1": 0.0},
    "E": {"e0": -3.63, "e1": 0.0},
    "A": {"e0": -1.76, "e1": 0.0},
    "T": {"e0": -1.47, "e1": 0.0},
    "Q": {"e0": -1.25, "e1": 0.00},
    "D": {"e0": -4.36, "e1": 0.0},
    "G": {"e0": -0.10, "e1": 0.00},
    "S": {"e0": -0.23, "e1": 0.00},
    "N": {"e0": -3.40, "e1": 0.0},
    "R": {"e0": -5.72, "e1": 0.0},
    "K": {"e0": -5.61, "e1": 0.0},
    "H": {"e0": -0.74, "e1": 0.0},
}

# Alkyl chains in a reversed-phase liquid chromatography stationary phase
# DeVido et al. 1998
energies_alkyl_chains = {
    "F": {"e0": -12.3, "e1": 0.0},
    "I": {"e0": -11.4, "e1": 0.0},
    "W": {"e0": -11.3, "e1": 0.0},
    "L": {"e0": -8.95, "e1": 0.0},
    "Y": {"e0": -8.70, "e1": 0.0},
    "M": {"e0": -8.23, "e1": 0.0},
    "V": {"e0": -7.80, "e1": 0.0},
    "P": {"e0": -7.15, "e1": 0.0},
    "C": {"e0": -4.05, "e1": 0.0},
    "E": {"e0": -3.0, "e1": 0.0},
    "A": {"e0": -2.67, "e1": 0.0},
    "T": {"e0": -2.50, "e1": 0.0},
    "Q": {"e0": -1.95, "e1": 0.0},
    "D": {"e0": -1.03, "e1": 0.0},
    "G": {"e0": -0.10, "e1": 0.00},
    "S": {"e0": -0.42, "e1": 0.00},
    "N": {"e0": -1.25, "e1": 0.00},
    "R": {"e0": -1.75, "e1": 0.00},
    "K": {"e0": -3.00, "e1": 0.00},
    "H": {"e0": -4.22, "e1": 0.00},
}


AA3_FROM_1 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def energies_ref(energies_ref="bulky"):
    ref_table = energies_alkyl_chains if energies_ref == "alkyl" else energies_bulk_octabol
    return {
        AA3_FROM_1[k]: {"e1": v["e1"], "e0": v["e0"]}
        for k, v in ref_table.items()
    }


# Kyte-doolite https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
KD = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9,
    "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
    "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

AA3 = set(KD.keys())

CHARGE = {aa: 0.0 for aa in AA3}
CHARGE.update({"ASP": -1.0, "GLU": -1.0, "LYS": +1.0, "ARG": +1.0, "HIS": +0.1})

VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
}
DEFAULT_VDW_RADIUS = 1.70


def _selected_residues(
    u: mda.Universe,
    protein_only: bool = True,
    chain: str | None = "A",
):
    residues = []
    for res in u.residues:
        chain_label = _residue_chain_label(res)
        if chain is not None and chain_label != chain:
            continue
        if protein_only and res.resname.strip() not in AA3:
            continue
        residues.append(res)
    return residues


def _infer_element(atom) -> str:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element

    atom_type = str(getattr(atom, "type", "") or "").strip().upper()
    if atom_type:
        return atom_type

    atom_name = "".join(ch for ch in str(atom.name).upper() if ch.isalpha())
    if atom_name.startswith("CL"):
        return "CL"
    if atom_name:
        return atom_name[0]
    return "C"


def _atom_radius(atom) -> float:
    return float(VDW_RADII.get(_infer_element(atom), DEFAULT_VDW_RADIUS))


def _residue_atom_data(residue) -> Tuple[np.ndarray, np.ndarray]:
    return residue.atoms.positions.astype(float), np.asarray(
        [_atom_radius(a) for a in residue.atoms], dtype=float
    )


def _residue_min_distance_and_contact(
    xyz_i: np.ndarray,
    r_i: np.ndarray,
    xyz_j: np.ndarray,
    r_j: np.ndarray,
    radius_buffer: float,
) -> Tuple[float, bool]:
    deltas = xyz_i[:, None, :] - xyz_j[None, :, :]
    dists = np.sqrt(np.sum(deltas * deltas, axis=2))
    min_dist = float(np.min(dists))
    contact_threshold = r_i[:, None] + r_j[None, :] + radius_buffer
    has_contact = bool(np.any(dists <= contact_threshold))
    return min_dist, has_contact


# MDAnalysis helper
def residue_coms(u: mda.Universe, protein_only: bool = True, chain="A") -> np.ndarray:
    coords = []
    for res in u.residues:
        if chain is not None and res.segid.strip() != chain:
            continue
        if protein_only and res.resname.strip() not in AA3:
            continue
        coords.append(res.atoms.center_of_mass())
    return np.array(coords, float)


def build_peptide_pairs(n_res: int) -> List[Tuple[int, int]]:
    return [(i, i + 1) for i in range(n_res - 1)]


def build_disulfide_pairs(
    u: mda.Universe,
    protein_only: bool = True,
    chain: str | None = "A",
    max_sg_distance: float = 2.2,
) -> List[Tuple[int, int]]:
    residues = _selected_residues(u, protein_only=protein_only, chain=chain)
    cys_sg = []
    for idx, res in enumerate(residues):
        if res.resname.strip() != "CYS":
            continue
        sg_atoms = res.atoms.select_atoms("name SG")
        if len(sg_atoms) == 0:
            continue
        cys_sg.append((idx, sg_atoms.positions[0].astype(float)))

    pairs: List[Tuple[int, int]] = []
    for i in range(len(cys_sg)):
        idx_i, pos_i = cys_sg[i]
        for j in range(i + 1, len(cys_sg)):
            idx_j, pos_j = cys_sg[j]
            if float(np.linalg.norm(pos_i - pos_j)) <= max_sg_distance:
                pairs.append((idx_i, idx_j))
    return pairs


def _contact_weight_35(
    xyz_i: np.ndarray,
    xyz_j: np.ndarray,
    box=None,
    eps: float = 1e-6,
) -> float:
    d = distance_array(xyz_i, xyz_j, box=box)
    d = d[d < 3.5]

    if d.size == 0:
        return 0.0

    return float(np.sum(1.0 / np.maximum(d, eps)))


def _residues_are_covalently_bonded(res_i, res_j) -> bool:
    """Check if any atom in res_i shares a bond with any atom in res_j."""
    atoms_j = set(res_j.atoms.indices)
    for atom in res_i.atoms:
        if not hasattr(atom, "bonded_atoms"):
            return False  # no bond information available
        for bonded in atom.bonded_atoms:
            if bonded.index in atoms_j:
                return True
    return False


# Legacy contact-weight pairs. The main pipeline now uses build_Kij_pairs.
def build_vdW_pairs(
    u: mda.Universe,
    protein_only: bool = True,
    chain: str | None = "A",
    exclude_bonded: bool = True,
    min_seq_separation: int = 1,
    use_pbc: bool = False,
) -> Tuple[List[Tuple[int, int]], Dict[Tuple[int, int], float]]:

    residues = _selected_residues(
        u,
        protein_only=protein_only,
        chain=chain,
    )

    atom_xyz = [
        res.atoms.positions.astype(np.float64)
        for res in residues
    ]

    box = u.dimensions if use_pbc else None

    pairs: List[Tuple[int, int]] = []
    weights: Dict[Tuple[int, int], float] = {}

    for i, xyz_i in enumerate(atom_xyz):
        for j in range(i + 1, len(atom_xyz)):
            if exclude_bonded and (j - i) <= min_seq_separation:
                continue

            w = _contact_weight_35(
                xyz_i,
                atom_xyz[j],
                box=box,
            )

            if w > 0:
                pair = (i, j)
                pairs.append(pair)
                weights[pair] = w

    return pairs, weights


# Excitation & energy model
def sigmoid(x, bias=-2.4):
    return 1.0 / (1.0 + np.exp(-(x + bias)))


def norm01(x: np.ndarray) -> np.ndarray:
    m = float(np.max(x)) if np.max(x) > 0 else 1.0
    return (x / m).astype(float)


# Legacy hydrophobicity scoring helper.
def kd_to_scores(resnames: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    kd_vals = np.array([KD.get(r, 0.0) for r in resnames], dtype=float)
    kd_min, kd_max = -4.5, 4.5
    kd_norm = (kd_vals - kd_min) / (kd_max - kd_min)   # [0,1] hydrophobic
    hydrophobic = np.clip(kd_norm, 0.0, 1.0)
    hydrophilic = 1.0 - hydrophobic
    return hydrophobic, hydrophilic


# Legacy hydrophobicity classifier.
def is_hydrophobic(resnames: List[str]):
    kd_vals = np.array([KD.get(r, 0.0) for r in resnames], dtype=float)
    return kd_vals > 0


# RY rotation (amplitude)
def angles_from_probability(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 0.0, 1.0)
    return 2.0 * np.arcsin(np.sqrt(p))


def _probability_to_ry_angle(
    p: Sequence[float] | float,
) -> np.ndarray:
    """
    Vectorized probability-to-RY angle conversion.
    """
    p = np.asarray(p, dtype=float)
    p = np.clip(p, 0.0, 1.0)

    return 2.0 * np.arcsin(np.sqrt(p))


def _residue_chain_label(residue) -> str:
    chain_label = str(getattr(residue, "chainID", "") or "").strip()
    if chain_label:
        return chain_label
    try:
        chain_ids = getattr(residue.atoms, "chainIDs", None)
        if chain_ids is not None and len(chain_ids) > 0:
            chain_label = str(chain_ids[0] or "").strip()
            if chain_label:
                return chain_label
    except Exception:
        pass
    return str(getattr(residue, "segid", "") or "").strip()


def boltzmann_expectation(delta_E, T=298.0, k_B=8.314462618e-3):
    """
    Compute expected energy, probabilities, and the probability
    at the expected energy level for a two-level Boltzmann system.
    """
    # Here the value is set to k_B = R; the gas constant
    beta = 1.0 / (k_B * T)
    Z = 1.0 + np.exp(-beta * delta_E)

    # State probabilities
    P0 = 1.0 / Z
    P1 = np.exp(-beta * delta_E) / Z

    # Energies
    E0, E1 = 0.0, delta_E

    # Expected energy
    E_expect = E0 * P0 + E1 * P1

    # Probability evaluated at expected energy
    P_at_Eexpect = np.exp(-beta * E_expect) / Z

    return {
        "E_expect": E_expect,
        "E0": E0,
        "E1": E1,
        "P0": P0,
        "P1": P1,
        "P(E_expect)": P_at_Eexpect,
    }


# protein propensities
@debug_mode
def excitation_propensities(SASA, etype=1):

    propensities = {}
    ref_energies = energies_ref("bulky")
    print("ref", ref_energies)
    for k, v in SASA.items():
        energy = ref_energies[k[2]]["e1"] - ref_energies[k[2]]["e0"]
        if KD[k[2]] > 0.5:
            propensities[k] = boltzmann_expectation(energy - v * energy)["P1"]
        elif KD[k[2]] < 0.5 and KD[k[2]] > -0.5:
            if v >= 0.1:
                propensities[k] = boltzmann_expectation(energy - v * energy)["P1"]
            else:
                v = 0.1
                propensities[k] = boltzmann_expectation(energy - v * energy)["P1"]
        else:
            propensities[k] = boltzmann_expectation(energy - (1 - v) * energy)["P1"]
    return propensities


def _make_sampler_and_backend(
    name: str,
    shots: int = 100_000,
    seed: int | None = None,
) -> Tuple[BaseSamplerV2, Any | None]:
    """
    Build a Qiskit v2.6.2 SamplerV2 primitive and optional BackendV2.

    Supported names:
      - "statevector_sampler", "default.qubit", "cpu":
            Pure Qiskit StatevectorSampler (no qiskit-aer dependency).
      - "aer_cpu", "aer_statevector":
            qiskit_aer.AerSimulator(method="statevector", device="CPU")
            wrapped in BackendSamplerV2.
      - "aer_gpu", "lightning.gpu", "gpu":
            qiskit_aer.AerSimulator(method="statevector", device="GPU")
            wrapped in BackendSamplerV2.
    """
    normalized = (name or "statevector_sampler").strip().lower()

    if normalized in {"statevector_sampler", "default.qubit", "cpu"}:
        return StatevectorSampler(default_shots=shots, seed=seed), None

    if normalized in {"aer_cpu", "aer_statevector", "aer_gpu", "lightning.gpu", "gpu"}:
        from qiskit_aer import AerSimulator

        device = "GPU" if normalized in {"aer_gpu", "lightning.gpu", "gpu"} else "CPU"
        aer_backend = AerSimulator(method="statevector", device=device)
        options: Dict[str, Any] = {"default_shots": shots}
        if seed is not None:
            options["seed_simulator"] = seed
        sampler = BackendSamplerV2(backend=aer_backend, options=options)
        return sampler, aer_backend

    # Mirror original make_device fallback for unrecognized backend names.
    return StatevectorSampler(default_shots=shots, seed=seed), None


def make_device(
    n_qubits: int,
    name: str,
    shots: int = 100_000,
    seed: int | None = None,
) -> BaseSamplerV2:
    """
    Backward-compatible device/sampler factory returning a Qiskit BaseSamplerV2.
    (n_qubits is accepted for signature parity with the PennyLane implementation.)
    """
    _ = n_qubits
    try:
        sampler, _ = _make_sampler_and_backend(name, shots=shots, seed=seed)
        return sampler
    except Exception:
        return StatevectorSampler(default_shots=shots, seed=seed)


def _lookup_kij(
    kijs: Mapping[Tuple[int, int], float],
    i: int,
    j: int,
) -> float | None:
    """
    Look up K_ij from a symmetric dictionary.

    Supports both:

        kijs[(i, j)]

    and

        kijs[(j, i)].
    """
    if (i, j) in kijs:
        return float(kijs[(i, j)])

    if (j, i) in kijs:
        return float(kijs[(j, i)])

    return None


def compute_directed_transfer_schedule(
    exc_propensities: Sequence[float],
    pairs: Sequence[Tuple[int, int]],
    kijs: Mapping[Tuple[int, int], float],
    E0: Sequence[float],
    E1: Sequence[float],
    beta: float = 1.0,  # \beta ~ 4.01 x 10^-4 mol J^-1 ; and \delta G is in kj mol-1;
    eps: float = 1.0e-12,
    gamma: float = 0.95,
    max_layers: int | None = None,
    min_deltaP: float = 0.0,
    phi_scale: float = 1.0,
    phi_cap: float | None = np.pi / 2,
    discount_start_power: int = 1,
) -> List[Dict[str, Any]]:
    """
    Compute the CRY entanglement schedule from directed residue-residue
    excitation transfer probabilities.
    """
    p = np.asarray(exc_propensities, dtype=float).reshape(-1)
    E0 = np.asarray(E0, dtype=float).reshape(-1)
    E1 = np.asarray(E1, dtype=float).reshape(-1)

    n_qubits = len(p)

    if len(E0) != n_qubits or len(E1) != n_qubits:
        raise ValueError(
            "E0, E1, and exc_propensities must have the same length."
        )

    if not (0.0 < gamma <= 1.0):
        raise ValueError("gamma must satisfy 0 < gamma <= 1.")

    p_clip = np.clip(p, eps, 1.0 - eps)

    dG = E1 - E0

    mu_values = np.array(
        [
            beta * float(dG[i]) * p_clip[i]
            for i in range(n_qubits)
        ],
        dtype=float,
    )

    # Build undirected neighbour list from pairs and K_ij values.
    neighbours: List[List[Tuple[int, float]]] = [
        [] for _ in range(n_qubits)
    ]

    for pair in pairs:
        i, j = int(pair[0]), int(pair[1])

        if i == j:
            continue

        if i < 0 or j < 0 or i >= n_qubits or j >= n_qubits:
            continue

        Kij = _lookup_kij(kijs, i, j)

        if Kij is None:
            continue

        Kij = float(Kij)

        if not np.isfinite(Kij) or Kij <= 0.0:
            continue

        neighbours[i].append((j, Kij))
        neighbours[j].append((i, Kij))

    # Denominator for each source residue i
    denom = np.zeros(n_qubits, dtype=float)

    for i in range(n_qubits):
        total = 0.0
        for k, Kik in neighbours[i]:
            total += Kik
        denom[i] = total

    # Compute all directed transfer probabilities DeltaP_{i -> j}.
    transfers: List[Dict[str, Any]] = []

    for i in range(n_qubits):
        if denom[i] <= 0.0:
            continue

        for j, Kij in neighbours[i]:
            drive_ij = max(
                0.0,
                float((mu_values[i] - mu_values[j]) / (mu_values[i] + mu_values[j] + 1e-5)),
            )

            if drive_ij <= 0.0:
                continue

            numerator = 2 * float(Kij) * drive_ij
            deltaP = numerator / (denom[i] + denom[j] + eps)
            deltaP = float(np.clip(deltaP, 0.0, 1.0))

            if deltaP <= min_deltaP:
                continue

            transfers.append(
                {
                    "i": int(i),
                    "j": int(j),
                    "Kij": float(Kij),
                    "mu_i": float(mu_values[i]),
                    "mu_j": float(mu_values[j]),
                    "drive": float(drive_ij),
                    "denom_i": float(denom[i]),
                    "deltaP": float(deltaP),
                }
            )

    # Highest directed transfer goes first.
    transfers.sort(
        key=lambda item: item["deltaP"],
        reverse=True,
    )

    if max_layers is not None:
        transfers = transfers[: int(max_layers)]

    # Assign layers and discounted CRY angles.
    for layer_index, entry in enumerate(transfers, start=1):
        deltaP = float(entry["deltaP"])
        power = (layer_index - 1) + int(discount_start_power)
        discount = float(gamma ** power)
        phi = float(_probability_to_ry_angle(discount * deltaP))

        if phi_cap is not None:
            phi = float(np.clip(phi, -abs(phi_cap), abs(phi_cap)))

        entry["layer"] = int(layer_index)
        entry["discount_power"] = int(power)
        entry["discount"] = float(discount)
        entry["phi"] = float(phi)

    return transfers


def _safe_logit_probability(
    x: float,
    eps: float = 1.0e-12,
) -> float:
    """
    Numerically stable logit:

        logit(x) = log[x / (1 - x)]

    with clipping to avoid log(0).
    """
    x = float(np.clip(x, eps, 1.0 - eps))
    return float(np.log(x / (1.0 - x)))


@debug_mode
def prepare_state_propensity(
    exc_propensities,
    kij_pairs,
    kijs,
    E0,
    E1,
    yy_cap: float | None = np.pi / 6.0,
    apply_yy_diffusion: bool = True,
    return_schedule: bool = True,
    yy_probability_threshold: float = 0.02,
    schedule: List[Dict[str, Any]] | None = None,
    qc: QuantumCircuit | None = None,
) -> QuantumCircuit | Tuple[QuantumCircuit, List[Dict[str, Any]]]:
    """
    Build (or append to) the thermodynamic excitation-redistribution QuantumCircuit
    driven by Kij couplings.

    Gate correspondence with PennyLane 0.45.1:
      - qml.RY(theta, wires=i)            -> qc.ry(theta, i)
      - qml.IsingYY(theta, wires=[i, j])  -> qc.ryy(theta, i, j)
      - qml.CRY(phi, wires=[i, j])        -> qc.cry(phi, i, j)

    Returns
    -------
    (qc, schedule) if return_schedule=True, else qc.
    """
    # 1) Local state initialization via RY rotations
    theta_local = angles_from_probability(np.asarray(exc_propensities, dtype=float))
    n_qubits = len(theta_local)
    if n_qubits == 0:
        raise ValueError("exc_propensities must be non-empty.")

    if qc is None:
        qc = QuantumCircuit(n_qubits)
    elif qc.num_qubits < n_qubits:
        raise ValueError(
            f"Provided QuantumCircuit has {qc.num_qubits} qubits, "
            f"but {n_qubits} are required."
        )

    for i in range(n_qubits):
        qc.ry(float(theta_local[i]), i)

    # 2) Classical computation of directed transfer layers
    if schedule is None:
        schedule = compute_directed_transfer_schedule(
            exc_propensities=exc_propensities,
            pairs=kij_pairs,
            kijs=kijs,
            E0=E0,
            E1=E1,
        )

    print("----------- Entaglement Layer -------\n")
    print(f"Schedule: {schedule}\n")
    print("------------  ")

    # 3) Apply CRY or RYY layers in descending DeltaP order
    for entry in schedule:
        i = int(entry["i"])
        j = int(entry["j"])
        phi = float(entry["phi"])

        # Use drive as the weak-diffusion proxy if available.
        deltaP_eff = float(
            entry.get(
                "drive",
                entry.get("deltaP", 0.0),
            )
        )

        # Weak diffusion correction: IsingYY(theta) == RYY(theta) = exp(-i * theta/2 * (Y x Y))
        if apply_yy_diffusion and deltaP_eff <= yy_probability_threshold:
            theta_yy = float(_probability_to_ry_angle(deltaP_eff))

            if yy_cap is not None:
                theta_yy = float(
                    np.clip(
                        theta_yy,
                        -abs(yy_cap),
                        abs(yy_cap),
                    )
                )

            if abs(theta_yy) > 0.0:
                qc.ryy(theta_yy, i, j)
        else:
            qc.cry(phi, i, j)
            qc.cry(-phi, j, i)

    if return_schedule:
        return qc, schedule

    return qc


def sample_bits(
    exc_propensities: np.ndarray,
    kij_pairs: Sequence[Tuple[int, int]],
    kijs: Dict[Tuple[int, int], float],
    E0,
    E1,
    cfg: Config,
    return_schedule: bool = False,
) -> np.ndarray | tuple[np.ndarray, List[Dict[str, Any]]]:
    n_qubits = len(exc_propensities)
    if n_qubits == 0:
        raise ValueError("exc_propensities must be non-empty.")

    schedule = compute_directed_transfer_schedule(
        exc_propensities=exc_propensities,
        pairs=kij_pairs,
        kijs=kijs,
        E0=E0,
        E1=E1,
    )

    qc = prepare_state_propensity(
        exc_propensities,
        kij_pairs,
        kijs,
        E0,
        E1,
        return_schedule=False,
        schedule=schedule,
    )
    meas = ClassicalRegister(n_qubits, "meas")
    qc.add_register(meas)
    qc.measure(range(n_qubits), meas)

    device_name = cfg.backend if cfg.backend is not None else cfg.pl_device

    def _run_sampling(target_device: str) -> np.ndarray:
        sampler, aer_backend = _make_sampler_and_backend(
            target_device,
            shots=cfg.shots,
            seed=cfg.seed,
        )
        if aer_backend is not None:
            pm = generate_preset_pass_manager(optimization_level=0, backend=aer_backend)
            qc_exec = pm.run(qc)
        else:
            qc_exec = qc

        job = sampler.run([qc_exec], shots=cfg.shots)
        pub_result = job.result()[0]

        # order="little" maps column i -> classical bit i -> qubit i,
        # matching PennyLane's [qml.sample(qml.PauliZ(i)) for i in range(n_qubits)].
        bits = pub_result.data.meas.to_bool_array(order="little")
        # Pauli-Z eigenvalues: |0> (False) -> +1, |1> (True) -> -1
        return np.where(bits, -1, 1)

    try:
        samples = _run_sampling(device_name)
    except Exception as exc:
        err_text = str(exc).lower()
        normalized = (device_name or "").strip().lower()
        should_fallback = normalized in {"lightning.gpu", "aer_gpu", "gpu"} and (
            "custatevec" in err_text
            or "lightning" in err_text
            or "cuda" in err_text
            or "gpu" in err_text
            or "qiskit_aer" in err_text
            or "aer" in err_text
        )
        if not should_fallback:
            raise
        print(
            f"Warning: {device_name} execution failed ({exc}). "
            "Falling back to statevector_sampler."
        )
        samples = _run_sampling("statevector_sampler")

    if return_schedule:
        return samples, schedule
    return samples


# --------------------- Energy accumulation ------------------
@debug_mode
def energy_from_samples(
    Z_samples: np.ndarray,
    E0: np.ndarray,
    E1: np.ndarray,
) -> np.ndarray:
    """
    Compute total configuration energy for each quantum sample.
    Ensures physical consistency with theoretical bounds.
    """
    ground = (Z_samples == +1)  # Indicator: Z-samples == +1
    # base energy per configuration
    E_res = ground * E0 + (~ground) * E1
    E_total = np.sum(E_res, axis=1)
    return E_total


def per_residue_energies(resnames):
    ref = energies_ref()
    return (
        np.array([ref[each]["e0"] for each in resnames]),
        np.array([ref[each]["e1"] for each in resnames]),
    )
```

---

## 4. Before-and-After Verification Test Suite

Save the original PennyLane file as `circuit_pl.py` and the ported Qiskit file as `circuit_qk.py` (or adjust the two imports at the top). The test fixture constructs a 4-residue system (`TRP`, `ILE`, `PHE`, `GLY`) that triggers **both** the paired `CRY` pathway ($\text{drive} > 0.02$) and the weak `IsingYY`/`RYY` diffusion pathway ($0 < \text{drive} \le 0.02$), verifying:

1. **Schedule identity** (exact match on every field).
2. **Exact statevector & full probability distribution identity** ($\text{atol} = 10^{-12}$, zero shot noise) using Qiskit's `Statevector.reverse_qargs()` to align little-endian Qiskit with big-endian PennyLane.
3. **Deterministic endianness & Pauli-Z sign check** (catches any bit-reversal or sign flip immediately).
4. **High-shot end-to-end sampling & `energy_from_samples` consistency** against both analytic expectation values and cross-framework results.
5. **Legacy `Config` backward compatibility** (`pl_device="default.qubit"` and `"lightning.gpu"`).

```python
"""
Before-and-After Equivalence Test: PennyLane 0.45.1 vs Qiskit v2.6.2
====================================================================
Run with:
    pytest -v test_port_equivalence.py
or directly:
    python test_port_equivalence.py
"""
import numpy as np
import pennylane as qml
from qiskit.quantum_info import Statevector

import circuit_pl  # Original PennyLane 0.45.1 circuit.py
import circuit_qk  # Ported Qiskit v2.6.2 circuit.py


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
```