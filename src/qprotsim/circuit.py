"""
Quantum Energy Propensity for Proteins — KE-aware & ground-skewed
=================================================================

- Qubit residue r: |0> = stable solvation, |1> = excited solvation.
- p_exc(r) = normalised(SASA_i) or 1-normalised(SASA_i) -> depending on hydrophobicity
- E0(r) negative; E1(r) = E0(r) + gap(r). The gap depends on the favourable solvation energy.
- Residue-residue Kij coupling factors drive directed excitation-transfer layers.
- Strong directed transfers use paired CRY gates; weak transfers can use YY diffusion.

Author: Pratik
"""
from dataclasses import dataclass, field
from numbers import Integral
import sys
from typing import Dict, List, Sequence, Tuple, Mapping, Any
import MDAnalysis as mda
import numpy as np
import pennylane as qml
from MDAnalysis.lib.distances import distance_array

from qprotsim.info_theory import debug_mode

@dataclass
class Config:
    """Structure and sampling settings, using PennyLane on the CPU by default.

    ``backend`` selects ``pennylane`` or native ``qiskit`` execution.
    ``backend_options`` are device options for the selected framework. Samples
    from either framework use Pauli-Z eigenvalues, never computational bits.
    """

    pdb_path: str
    protein_only: bool = True

    # shots and device
    shots: int = 100_000
    pl_device: str = "default.qubit"
    backend: str = "pennylane"
    qiskit_backend: str = "aer_simulator"
    backend_options: dict[str, Any] = field(default_factory=dict)
    seed: int | None = None

    def __post_init__(self):
        if self.backend not in ("pennylane", "qiskit"):
            raise ValueError("backend must be 'pennylane' or 'qiskit'.")
        if isinstance(self.shots, bool) or not isinstance(self.shots, Integral) or self.shots <= 0:
            raise ValueError("shots must be a positive integer.")
        if not isinstance(self.pl_device, str) or not self.pl_device.strip():
            raise ValueError("pl_device must be a non-empty PennyLane device name.")
        if not isinstance(self.qiskit_backend, str) or not self.qiskit_backend.strip():
            raise ValueError("qiskit_backend must be a non-empty Qiskit Aer backend name.")
        if not isinstance(self.backend_options, Mapping) or not all(
            isinstance(key, str) for key in self.backend_options
        ):
            raise ValueError("backend_options must be a mapping with string keys.")
        self.backend_options = dict(self.backend_options)
        if self.seed is not None and (
            isinstance(self.seed, bool) or not isinstance(self.seed, Integral) or self.seed < 0
        ):
            raise ValueError("seed must be a non-negative integer or null.")

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
    "G": {"e0":  -0.10, "e1":  0.00},
    "S": {"e0": -0.23, "e1": 0.00},
    "N": {"e0": -3.40, "e1": 0.0},
    "R": {"e0": -5.72, "e1": 0.0},
    "K": {"e0": -5.61, "e1": 0.0},
    "H": {"e0": -0.74, "e1": 0.0}
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
    "G": {"e0":  -0.10, "e1": 0.00},
    "S": {"e0": -0.42, "e1": 0.00},
    "N": {"e0": -1.25, "e1": 0.00},
    "R": {"e0": -1.75, "e1": 0.00},
    "K": {"e0": -3.00, "e1": 0.00},
    "H": {"e0": -4.22, "e1": 0.00}
}


AA3_FROM_1 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
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
    "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5
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
    return residue.atoms.positions.astype(float), np.asarray([_atom_radius(a) for a in residue.atoms], dtype=float)


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
def sigmoid(x, bias = -2.4): return 1.0 / (1.0 + np.exp(-(x + bias)))

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
def is_hydrophobic(resnames:List[str]):
    kd_vals = np.array([KD.get(r, 0.0) for r in resnames], dtype=float)
    return kd_vals > 0

# RY rotation (amplitude)
def angles_from_probability(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 0.0, 1.0)
    return 2.0*np.arcsin(np.sqrt(p))

def _probability_to_ry_angle(
    p: Sequence[float],
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
        "P(E_expect)": P_at_Eexpect
    }

# protein propensities
@debug_mode
def excitation_propensities(SASA, etype = 1):

    propensities = {}
    ref_energies = energies_ref("bulky")
    print("ref",ref_energies)
    for k,v in SASA.items():
        energy = ref_energies[k[2]]["e1"] - ref_energies[k[2]]["e0"]
        if KD[k[2]] > 0.5:
            propensities[k] = boltzmann_expectation( energy - v*energy)["P1"]
        elif KD[k[2]] < 0.5 and KD[k[2]] > -0.5:
            if v >= 0.1:
                propensities[k] = boltzmann_expectation( energy - v*energy)["P1"]
            else:
                v = 0.1
                propensities[k] = boltzmann_expectation(energy - v*energy)["P1"]
        else:
            propensities[k] = boltzmann_expectation(energy - (1-v)*energy)["P1"]
    return propensities

def make_device(n_qubits: int, name: str, **device_options):
    try:
        return qml.device(name, wires=n_qubits, **device_options)
    except Exception:
        # Explicit options must fail visibly rather than silently being ignored.
        if device_options:
            raise
        return qml.device("default.qubit", wires=n_qubits)

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
    beta: float = 1.0, # \beta ~ 4.01 x 10^-4 mol J^-1 ; and \delta G is in kj mol-1;
    eps: float = 1.0e-12,
    gamma: float = 0.95,
    max_layers: int | None = None,
    min_deltaP: float = 0.0,
    phi_scale: float = 1.0,
    phi_cap: float | None = np.pi/2,
    discount_start_power: int = 1,
) -> List[Dict[str, Any]]:
    """
    Compute the CRY entanglement schedule from directed residue-residue
    excitation transfer probabilities.

    For the clipped initial excitation probability P_i, the implementation uses

        mu_i = beta * (E1_i - E0_i) * P_i
        degree_i = sum_{k in N(i)} K_ik
        drive_ij = max(0, (mu_i - mu_j) / (mu_i + mu_j + 1e-5))
        DeltaP_{i -> j} = clip(
            2 * K_ij * drive_ij / (degree_i + degree_j + eps), 0, 1
        )

    DeltaP is a schedule score. Marginal excitation changes are determined by
    the full ordered circuit.

    The CRY angle for layer ell is:

        deltaP_discounted = gamma^ell * deltaP_raw
        phi = 2 arcsin(sqrt(deltaP_discounted))

    Parameters
    ----------
    exc_propensities
        Initial excitation probabilities P_i for each residue/qubit.

    pairs
        Undirected residue-residue pairs.

    kijs
        Dictionary mapping residue pairs to K_ij values.

    E0, E1
        Residue-wise solvation/free-energy values. DeltaG = E1 - E0.

    beta
        Numerical energy scale in mu_i, defaulting to 1.0.

    eps
        Small numerical stabilizer.

    gamma
        Layer discount factor. If gamma=0.99, then later layers receive
        factors gamma^1, gamma^2, gamma^3, ...

    max_layers
        Optional maximum number of directed transfer layers to apply.

    min_deltaP
        Ignore transfers with DeltaP <= min_deltaP.

    phi_scale
        Reserved argument. The current implementation does not apply it.

    phi_cap
        Optional maximum absolute CRY angle. Use None for no cap.

    discount_start_power
        If 1, layer 1 uses gamma^1.
        If 0, layer 1 is undiscounted and uses gamma^0.

    Returns
    -------
    schedule
        Sorted list of dictionaries. Each entry contains:

            i, j, Kij, mu_i, mu_j, drive, denom_i,
            deltaP, layer, discount_power, discount, phi
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
            # beta * float(dG[i]) + _safe_logit_probability(p_clip[i], eps=eps)
            beta*float(dG[i])*p_clip[i]
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

    # ------------------------------------------------------------
    # Denominator for each source residue i:
    #
    # denom_i = sum_k K_ik [mu_i - mu_k]_+
    # ------------------------------------------------------------
    denom = np.zeros(n_qubits, dtype=float)

    for i in range(n_qubits):
        total = 0.0
        for k, Kik in neighbours[i]:
            # drive_ik = max(0.0, float((mu_values[i] - mu_values[k])/(mu)))
            total += Kik
            # total += float(Kik) * drive_ik

        denom[i] = total

    # Compute all directed transfer probabilities DeltaP_{i -> j}.
    transfers: List[Dict[str, Any]] = []

    for i in range(n_qubits):

        if denom[i] <= 0.0:
            continue

        for j, Kij in neighbours[i]:

            drive_ij = max(0.0, float((mu_values[i] - mu_values[j]) / (mu_values[i] + mu_values[j] + 1e-5)))

            if drive_ij <= 0.0:
                continue

            numerator = 2*float(Kij) * drive_ij

            # deltaP = (
            #     float(p_clip[i])
            #     * (1.0 - float(p_clip[j]))
            #     * numerator
            #     / (float(denom[i]) + eps)
            # )
            deltaP = (
                # (1.0 - float(p_clip[j]))
               numerator
                / (denom[i] + denom[j] + eps)
            )

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
        phi = _probability_to_ry_angle(discount*deltaP)
        

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

def _propensity_gate_sequence(
    exc_propensities,
    schedule,
    yy_cap=np.pi / 6.0,
    apply_yy_diffusion=True,
    yy_probability_threshold=0.02,
):
    """Yield the same ordered RY, CRY and YY rotations for both frameworks."""
    probabilities = np.asarray(exc_propensities, dtype=float)
    if probabilities.ndim != 1 or probabilities.size == 0:
        raise ValueError("exc_propensities must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(probabilities)):
        raise ValueError("exc_propensities must contain finite probabilities.")
    theta_local = angles_from_probability(probabilities)
    for i, theta in enumerate(theta_local):
        yield "ry", float(theta), (i,)

    for entry in schedule:
        i, j = int(entry["i"]), int(entry["j"])
        if i == j or not (0 <= i < len(theta_local) and 0 <= j < len(theta_local)):
            raise ValueError("Schedule entries must name two distinct, valid qubit indices.")
        phi = float(entry["phi"])
        deltaP_eff = float(entry.get("drive", entry.get("deltaP", 0.0)))
        if not np.isfinite(phi) or not np.isfinite(deltaP_eff):
            raise ValueError("Schedule angles and transfer probabilities must be finite.")
        if apply_yy_diffusion and deltaP_eff <= yy_probability_threshold:
            theta_yy = float(_probability_to_ry_angle(deltaP_eff))
            if yy_cap is not None:
                theta_yy = float(np.clip(theta_yy, -abs(yy_cap), abs(yy_cap)))
            if abs(theta_yy) > 0.0:
                yield "ryy", theta_yy, (i, j)
        else:
            yield "cry", phi, (i, j)
            yield "cry", -phi, (j, i)


@debug_mode
def prepare_state_propensity(
    exc_propensities,
    kij_pairs,
    kijs,
    E0, E1,
    yy_cap: float | None = np.pi / 6.0,
    apply_yy_diffusion=True, return_schedule: bool = True, yy_probability_threshold=0.02,
    schedule: List[Dict[str, Any]] | None = None,
):
    """Prepare the PennyLane propensity circuit with shared ordered gate layers.

    Residue i is qubit i: |0> is the ground solvation state (Z=+1),
    and |1> is the excited solvation state (Z=-1).
    """
    if schedule is None:
        schedule = compute_directed_transfer_schedule(
            exc_propensities=exc_propensities,
            pairs=kij_pairs,
            kijs=kijs,
            E0=E0,
            E1=E1,
        )
    print("----------- Entanglement Layer -------\n")
    print(f"Schedule: {schedule}\n")
    print("------------  ")
    for gate, angle, wires in _propensity_gate_sequence(
        exc_propensities, schedule, yy_cap, apply_yy_diffusion, yy_probability_threshold
    ):
        if gate == "ry":
            qml.RY(angle, wires=wires[0])
        elif gate == "cry":
            qml.CRY(angle, wires=wires)
        else:
            qml.IsingYY(angle, wires=wires)
    return schedule if return_schedule else None


def build_qiskit_circuit(
    exc_propensities,
    kij_pairs,
    kijs,
    E0, E1,
    yy_cap: float | None = np.pi / 6.0,
    apply_yy_diffusion=True,
    yy_probability_threshold=0.02,
    schedule: List[Dict[str, Any]] | None = None,
    measure: bool = True,
):
    """Build a native Qiskit circuit matching the PennyLane propensity circuit.

    Residue i maps to qubit i and, when measured, classical bit i. Qiskit
    displays kets and memory strings as |q_(n-1) ... q_0>; this display order
    does not change the residue mapping. |0> has Z=+1 and |1> has Z=-1.
    Set ``measure=False`` for statevector inspection, drawing, or further gates.
    Qiskit is optional; install QProtSim with the ``qiskit`` extra to use it.
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError as exc:
        raise ImportError(
            'Qiskit support requires the optional packages: pip install ".[qiskit]" '
            '(or pip install "qprotsim[qiskit]").'
        ) from exc

    probabilities = np.asarray(exc_propensities, dtype=float)
    if probabilities.ndim != 1 or probabilities.size == 0:
        raise ValueError("exc_propensities must be a non-empty one-dimensional array.")
    n_qubits = len(probabilities)
    if schedule is None:
        schedule = compute_directed_transfer_schedule(
            probabilities, kij_pairs, kijs, E0, E1
        )
    circuit = QuantumCircuit(n_qubits, n_qubits if measure else 0, name="qprotsim")
    circuit.metadata = {
        "sample_convention": "Pauli-Z: |0> = +1 (ground), |1> = -1 (excited)",
        "residue_qubit_mapping": "residue i -> qubit i -> classical bit i",
        "memory_order": "q_(n-1) ... q_0; reverse before forming residue columns",
    }
    for gate, angle, wires in _propensity_gate_sequence(
        probabilities, schedule, yy_cap, apply_yy_diffusion, yy_probability_threshold
    ):
        getattr(circuit, gate)(angle, *wires)
    if measure:
        circuit.measure(range(n_qubits), range(n_qubits))
    return circuit


def bits_to_z_samples(bits) -> np.ndarray:
    """Convert a (shots, residues) 0/1 array to signed Pauli-Z samples.

    The conversion is explicit: computational |0> maps to +1 (ground),
    computational |1> maps to -1 (excited), using Z = 1 - 2 * bit.
    An all-ones input means all excited; no encoding is guessed from values.
    """
    bits = np.asarray(bits)
    if bits.ndim != 2 or not np.all((bits == 0) | (bits == 1)):
        raise ValueError("Computational-basis samples must be a 2D array containing only 0 and 1.")
    return 1 - 2 * bits.astype(np.int8)


def qiskit_memory_to_z_samples(memory: Sequence[str], n_qubits: int) -> np.ndarray:
    """Decode shot memory in Qiskit's n-1 ... 0 order into residue columns.

    For example, memory '001' means residue/qubit 0 is excited, so its
    Pauli-Z row is [-1, +1, +1]. Register-separating whitespace is ignored.
    """
    if isinstance(n_qubits, bool) or not isinstance(n_qubits, Integral) or n_qubits <= 0:
        raise ValueError("n_qubits must be a positive integer.")
    rows = []
    for shot in memory:
        if not isinstance(shot, str):
            raise ValueError("Qiskit memory must contain binary strings.")
        bitstring = "".join(shot.split())
        if len(bitstring) != n_qubits or any(bit not in "01" for bit in bitstring):
            raise ValueError(f"Each Qiskit memory string must contain exactly {n_qubits} binary digits.")
        rows.append([int(bit) for bit in reversed(bitstring)])
    bits = np.asarray(rows, dtype=np.int8).reshape(-1, n_qubits)
    return bits_to_z_samples(bits)


def _sample_qiskit(exc_propensities, kij_pairs, kijs, E0, E1, cfg, schedule):
    try:
        from qiskit import transpile
        from qiskit_aer import Aer
    except ImportError as exc:
        raise ImportError(
            'Qiskit execution requires Qiskit Aer: pip install ".[qiskit]" '
            '(or pip install "qprotsim[qiskit]").'
        ) from exc
    circuit = build_qiskit_circuit(
        exc_propensities, kij_pairs, kijs, E0, E1, schedule=schedule
    )
    try:
        backend = Aer.get_backend(cfg.qiskit_backend)
        backend.set_options(**cfg.backend_options)
    except Exception as exc:
        raise ValueError(
            f"Cannot configure Qiskit backend {cfg.qiskit_backend!r}: {exc}"
        ) from exc
    seed = None if cfg.seed is None else int(cfg.seed)
    compiled = transpile(circuit, backend, seed_transpiler=seed)
    run_options = {"shots": int(cfg.shots), "memory": True}
    if seed is not None:
        run_options["seed_simulator"] = seed
    result = backend.run(compiled, **run_options).result()
    if not result.success:
        raise RuntimeError(f"Qiskit backend {cfg.qiskit_backend!r} execution failed: {result.status}")
    samples = qiskit_memory_to_z_samples(result.get_memory(compiled), len(exc_propensities))
    if samples.shape != (cfg.shots, len(exc_propensities)):
        raise RuntimeError("Qiskit returned a sample shape different from (shots, residues).")
    return samples


def sample_bits(exc_propensities: np.ndarray,
                kij_pairs: Sequence[Tuple[int, int]],
                kijs: Dict[Tuple[int, int], float],
                E0, E1,
                cfg: Config,
                return_schedule: bool = False,
                ) -> np.ndarray | tuple[np.ndarray, List[Dict[str, Any]]]:
    """Sample either framework and return (shots, residues) Pauli-Z values.

    |0> is +1 (ground), |1> is -1 (excited). Qiskit's measured bits are
    explicitly converted and reordered before entering the analysis pipeline.
    """
    n_qubits = len(exc_propensities)
    if n_qubits == 0:
        raise ValueError("At least one residue/qubit is required for sampling.")
    schedule = compute_directed_transfer_schedule(
        exc_propensities=exc_propensities,
        pairs=kij_pairs,
        kijs=kijs,
        E0=E0,
        E1=E1,
    )
    if cfg.backend == "qiskit":
        samples = _sample_qiskit(
            exc_propensities, kij_pairs, kijs, E0, E1, cfg, schedule
        )
        return (samples, schedule) if return_schedule else samples

    def _run_sampling(device_name: str):
        device_options = dict(cfg.backend_options)
        if cfg.seed is not None:
            device_options["seed"] = int(cfg.seed)
        dev = make_device(n_qubits, device_name, **device_options)

        @qml.set_shots(cfg.shots)
        @qml.qnode(dev)
        def qnode():
            prepare_state_propensity(
                exc_propensities,
                kij_pairs,
                kijs,
                E0,
                E1,
                schedule=schedule,
            )
            return [qml.sample(qml.PauliZ(i)) for i in range(n_qubits)]

        out = qnode()
        # Single-shot measurements can be scalars; preserve both axes.
        return np.asarray(out, dtype=np.int8).reshape(n_qubits, cfg.shots).T

    try:
        samples = _run_sampling(cfg.pl_device)
    except Exception as exc:
        err_text = str(exc).lower()
        should_fallback = (
            cfg.pl_device == "lightning.gpu"
            and (
                "custatevec" in err_text
                or "lightning" in err_text
                or "cuda" in err_text
                or "gpu" in err_text
            )
        )
        if not should_fallback:
            raise
        print(
            f"Warning: {cfg.pl_device} execution failed ({exc}). "
            "Falling back to default.qubit."
        )
        samples = _run_sampling("default.qubit")

    if return_schedule:
        return samples, schedule
    return samples


# --------------------- Energy accumulation ------------------
@debug_mode
def energy_from_samples(Z_samples: np.ndarray,
                        E0: np.ndarray,
                        E1: np.ndarray,
    ) -> np.ndarray:
    """
    Compute total energy from explicit Pauli-Z samples, shaped (shots, residues).

    |0> is +1 and contributes E0; |1> is -1 and contributes E1.
    Convert computational 0/1 bits with ``bits_to_z_samples`` first. Input
    values are never used to guess the convention: all +1 means all ground.
    """
    Z_samples = np.asarray(Z_samples)
    if Z_samples.ndim != 2 or not np.all((Z_samples == 1) | (Z_samples == -1)):
        raise ValueError(
            "Z_samples must be a 2D array of Pauli-Z values +1/-1. "
            "Convert 0/1 bits with bits_to_z_samples first."
        )
    E0, E1 = np.asarray(E0, dtype=float), np.asarray(E1, dtype=float)
    if E0.shape != (Z_samples.shape[1],) or E1.shape != E0.shape:
        raise ValueError("E0 and E1 must be one-dimensional with one energy per residue column.")
    ground = (Z_samples == +1) # Indicator: Z-samples == +1
    # print("Ground samples",ground,type(ground),ground.shape)
    # base energy per configuration
    E_res = ground * E0 + (~ground) * E1
    E_total = np.sum(E_res, axis=1)
    return E_total

def per_residue_energies(resnames):
    ref = energies_ref()
    return np.array([ref[each]["e0"] for each in resnames]), np.array([ref[each]["e1"] for each in resnames])
